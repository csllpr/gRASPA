# Performance Notes

This note documents the recent runtime changes that were added to improve throughput on modern NVIDIA GPUs without changing the default simulation path.

## Current opt-in path

The input keyword

```text
UseGPUReduction yes
```

now enables two GPU-side final-reduction paths:

1. Widom per-trial energy reduction
2. Ewald move-delta reduction for:
   `GPU_EwaldDifference_General`
   `GPU_EwaldDifference_IdentitySwap`
   `GPU_EwaldDifference_LambdaChange`

The default is still:

```text
UseGPUReduction no
```

That keeps the older host-side final summation behavior unless the user explicitly opts in.

## What changed internally

The recent runtime work in `src_clean/` did four things:

1. Moved blocked-pocket hot-path checks onto the GPU and replaced the old host scan with a small device summary.
2. Added device-side Widom trial-energy reduction so the host no longer needs the full per-block partial array for each Widom trial.
3. Added device-side Ewald delta reduction so the host copies back only the final same-type / cross-type pair when GPU reduction is enabled.
4. Reserved explicit tail scratch inside `Sim.Blocksum` and preserved that tail during `Blocksum` growth.

## Intended usage

Use `UseGPUReduction yes` when:

* the workload is Widom-heavy
* the workload uses charged systems with frequent Ewald move-delta evaluations
* many concurrent jobs are being packed onto one GPU through MPS

Leave it at `no` when:

* you want the historical default path
* you are comparing against older local runs and want to minimize runtime-path changes

## Validation status

The current implementation was checked against the designated example suite and the existing `pytest` example test.

For the preserved March 26, 2026 baseline runs, the cycle lines remained identical after the reduction changes.

## Benchmark observations

On the local GPU 1 benchmark setup used during this refactor work:

* the Widom GPU-reduction path gave the larger gain
* the added Ewald move-delta reduction was smaller but still positive
* removing explicit host/device rendezvous points in the hot path was not a win under many-process MPS load and was reverted

So the current recommendation is:

* keep explicit sync points where the host immediately consumes move results
* use `UseGPUReduction yes` for throughput-oriented Widom / charged workloads
* benchmark with your real workload before making it the default in production inputs

## Known limits

These changes do not solve the larger architectural limits by themselves:

* one single Monte Carlo chain still has a serial acceptance frontier
* the code still has remaining managed-memory / host-orchestration overhead outside these reduction paths
* external job concurrency can still saturate at a lower-throughput point if too many jobs share one GPU

The next likely optimization areas are:

1. further managed-memory cleanup
2. reducing remaining host-side move orchestration costs
3. larger-scale multi-walker or multi-chain execution inside one process

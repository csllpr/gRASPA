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
3. GPU-side final reduction for single-body / lambda-change VDW+real move energies

The default is still:

```text
UseGPUReduction no
```

That keeps the older host-side final summation behavior unless the user explicitly opts in.

There is also an opt-in host RNG path:

```text
UseFastHostRNG yes
```

The default remains:

```text
UseFastHostRNG no
```

That keeps the historical `std::rand()`-based host RNG behavior unless the user explicitly opts in to the faster generator.

## What changed internally

The recent runtime work in `src_clean/` did four things:

1. Moved blocked-pocket hot-path checks onto the GPU and replaced the old host scan with a small device summary.
2. Added device-side Widom trial-energy reduction so the host no longer needs the full per-block partial array for each Widom trial.
3. Added device-side Ewald delta reduction so the host copies back only the final same-type / cross-type pair when GPU reduction is enabled.
4. Added device-side final reduction for the single-body / CBCF lambda VDW+real move-energy path.
5. Reserved explicit tail scratch inside `Sim.Blocksum` and preserved that tail during `Blocksum` growth.
6. Replaced the DNN adsorbate-mask managed allocation with an explicit host pointer plus device mirror.
7. Replaced the `Vars.Sims` managed wrapper with a plain host allocation and removed the last kernels that depended on passing `Simulations` itself into device code.
8. Reduced small but repeated Widom-side CPU overhead by reusing the existing shifted-Boltzmann scratch vector, shrinking the first-bead trial-position kernel arguments, and removing the stale random-setup debug launch.
9. Added an opt-in faster host-side uniform RNG backend for CPU-side move selection and host random-buffer generation.
10. Packed the Widom overlap validity bit into the GPU-reduced per-trial payload so the `UseGPUReduction yes` path no longer needs a separate host copy of the trial flags.
11. Changed the opt-in Ewald move-delta path so the Fourier kernel can accumulate same-type and cross-type totals directly into the reduced scratch, removing the extra GPU reduction launch from `UseGPUReduction yes`.
12. Changed the opt-in Widom VDW+real trial-energy path so the main pair-energy kernel accumulates directly into the per-trial reduced payload, removing the extra GPU reduction launch from `UseGPUReduction yes`.
13. Increased the host/device random-buffer size used by `UseFastHostRNG yes` so the throughput-oriented path refills and copies the host random buffer less often.

## Intended usage

Use `UseGPUReduction yes` when:

* the workload is Widom-heavy
* the workload uses charged systems with frequent Ewald move-delta evaluations
* many concurrent jobs are being packed onto one GPU through MPS

Leave it at `no` when:

* you want the historical default path
* you are comparing against older local runs and want to minimize runtime-path changes

Use `UseFastHostRNG yes` when:

* the workload is CPU-limited by host-side move orchestration or random-buffer refill
* you are optimizing throughput rather than same-seed trajectory identity

Leave it at `no` when:

* you want the historical host RNG path
* you need to keep old same-seed trajectories comparable

## Validation status

The current implementation was checked against the designated example suite and the existing `pytest` example test.

For the preserved March 26, 2026 baseline runs, the cycle lines remained identical after the reduction changes.

## Benchmark observations

On the local GPU 1 benchmark setup used during this refactor work:

* the Widom GPU-reduction path gave the larger gain
* the added Ewald move-delta reduction was smaller but still positive
* the added single-body / lambda move-energy reduction was positive on a charged, move-heavy Bae benchmark
* `UseFastHostRNG yes` improved the standard 16-way Widom benchmark by about `1.0249x`
* packing the Widom validity bit into the reduced trial payload improved the 16-way `UseGPUReduction yes` + `UseFastHostRNG yes` benchmark by about `1.2052x` versus the immediately previous build
* direct Ewald pair accumulation in the opt-in path improved the same 16-way `UseGPUReduction yes` + `UseFastHostRNG yes` benchmark by about `1.0665x` versus the previous implementation
* direct Widom per-trial accumulation in the opt-in path improved the same 16-way `UseGPUReduction yes` + `UseFastHostRNG yes` benchmark by about `1.0492x` versus the previous implementation
* increasing the fast-RNG buffer size improved the same 16-way `UseGPUReduction yes` + `UseFastHostRNG yes` benchmark by about `1.0175x` versus the immediately previous build
* removing explicit host/device rendezvous points in the hot path was not a win under many-process MPS load and was reverted
* converting `Vars.Sims` to a host wrapper was throughput-neutral to slightly positive on the standard 16-way Widom benchmark while also removing the last `cudaMallocManaged(...)` dependency in `src_clean/`
* the Widom-side CPU-overhead cleanup produced a further small gain on the same 16-way benchmark

An `nsys` short-run profile on the throughput-oriented Widom benchmark also showed why this helped: the dominant CUDA API costs were call count, not payload size, with `cudaMemcpy`, `cudaDeviceSynchronize`, and `cudaLaunchKernel` accounting for almost all traced host runtime. Cutting one hot-path memcpy in the opt-in Widom reduction path therefore translated into a large throughput gain.

So the current recommendation is:

* keep explicit sync points where the host immediately consumes move results
* use `UseGPUReduction yes` for throughput-oriented Widom / charged workloads
* benchmark `UseFastHostRNG yes` on throughput-focused runs where same-seed identity is not required
* prefer the current `UseGPUReduction yes` path on Widom-heavy runs, since it now avoids one of the hottest per-move host copies
* benchmark with your real workload before making it the default in production inputs

## Known limits

These changes do not solve the larger architectural limits by themselves:

* one single Monte Carlo chain still has a serial acceptance frontier
* the code still has substantial host-orchestration overhead even after the managed-memory cleanup
* external job concurrency can still saturate at a lower-throughput point if too many jobs share one GPU

The next likely optimization areas are:

1. further managed-memory cleanup
2. reducing remaining host-side move orchestration costs
3. larger-scale multi-walker or multi-chain execution inside one process

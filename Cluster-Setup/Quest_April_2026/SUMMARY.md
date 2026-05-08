# Quest Installation Package Summary — April 2026

This directory contains a complete, ready-to-use installation package for compiling gRASPA (vanilla version) on Quest cluster at Northwestern University. The resulting binary runs on **both A100 and H100** GPU nodes.

## Package Contents

### Core Scripts
- **`NVC_COMPILE_QUEST_VANILLA`** — Main compilation script (A100 + H100)
- **`compile_graspa.job`** — SLURM job submission template

### Documentation
- **`README.md`** — Comprehensive installation guide with troubleshooting
- **`QUICKSTART.md`** — Quick start guide for experienced users
- **`FILES.md`** — File structure and overview
- **`INSTALLATION_CHECKLIST.md`** — Step-by-step checklist
- **`SUMMARY.md`** — This file

## Quick Overview

### What This Package Does
Compiles gRASPA vanilla version (without ML support) on Quest GPU nodes using the NVIDIA HPC SDK compiler, targeting both A100 (cc80) and H100 (cc90) GPUs.

### What You Need
- Quest account with GPU access
- gRASPA source code in `src_clean/` directory
- ~5 minutes for compilation

### What You Get
- Compiled executable: `nvc_main.x` (~5–6 MB)
- Works on both A100 and H100 nodes — use `--gres=gpu:1` (no need to specify type)

## Tested Configuration

| Component | Value |
|-----------|-------|
| **Cluster** | Quest, Northwestern University |
| **OS** | RHEL 8.10, kernel 4.18.0-553 |
| **NVIDIA Driver** | 570.86.15 (CUDA 12.8 capable) |
| **Compiler** | nvhpc/23.3-gcc-10.4.0 (spack_v20d1) |
| **CUDA Toolkit** | 11.8 (bundled in nvhpc 23.3) |
| **GPU Targets** | cc80 (A100) + cc90 (H100) |
| **A100 test** | qgpu2002 — exit code 0 |
| **H100 test** | qgpu3006 — exit code 0 |
| **Test case** | CO2 in Mg-MOF-74 (UFF), 10k init + 10k prod cycles |
| **Test date** | April 3, 2026 |

## Changes from Previous Quest Installation (Quest/)

| | Old (Quest/) | New (Quest_April_2026/) |
|--|-------------|------------------------|
| **Compiler** | nvhpc/21.9-gcc | nvhpc/23.3-gcc-10.4.0 |
| **CUDA toolkit** | 11.4 (default) | 11.8 (explicit `-gpu=cuda11.8`) |
| **GPU targets** | cc80 only (implicit) | cc80 + cc90 (explicit) |
| **H100 support** | No | Yes |
| **libatomic fix** | Not needed | gcc-10.3.0 lib64 prepended to linker path |
| **Runtime fix** | Not needed | Save/restore `CUDA_VISIBLE_DEVICES` around `module purge` |

## Known Issues on Quest (April 2026)

1. **gcc-10.4.0 libatomic.so** has compressed debug sections that system `ld 2.30` cannot read — the compile script works around this by prepending gcc-10.3.0's lib64 directory
2. **`module purge` clears `CUDA_VISIBLE_DEVICES`** — run scripts must save and restore this variable
3. **nvhpc/22.7-gcc** has an incompatible Thrust version — do not use
4. **nvhpc ≥24.x** also has the libatomic issue — use 23.3

## Next Steps After Installation

1. **Test the executable**: Run a gRASPA example from `Examples/`
2. **Create run scripts**: Use the template in README.md
3. **Configure resources**: Adjust memory/time limits based on your simulation needs

## Support

- Check `README.md` for detailed troubleshooting
- Review `INSTALLATION_CHECKLIST.md` for step-by-step verification
- Quest IT: https://www.it.northwestern.edu/departments/it-services-support/research/computing/quest/
- gRASPA docs: https://zhaoli2042.github.io/gRASPA-mkdoc

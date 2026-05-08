# gRASPA Installation on Quest — April 2026

This directory contains scripts and instructions for compiling gRASPA (vanilla version, without ML potential support) on the Quest cluster at Northwestern University. The resulting binary runs on **both A100 and H100** GPU nodes.

## Overview

gRASPA is a GPU-accelerated Monte Carlo simulation software for molecular adsorption in nanoporous materials. This installation guide covers the **vanilla version** without machine learning potential support.

## Prerequisites

- Access to Quest cluster with GPU nodes (A100 or H100)
- Quest account/project ID for job submission

## Quick Start

### Step 1: Prepare Source Code

Navigate to the gRASPA source directory:
```bash
cd /path/to/gRASPA/src_clean
```

### Step 2: Copy Compilation Script

Copy the compilation script to your source directory:
```bash
cp Cluster-Setup/Quest_April_2026/NVC_COMPILE_QUEST_VANILLA .
chmod +x NVC_COMPILE_QUEST_VANILLA
```

### Step 3: Compile gRASPA

You have two options:

#### Option A: Submit as a Job (Recommended)

1. Copy and edit the job submission script:
   ```bash
   cp Cluster-Setup/Quest_April_2026/compile_graspa.job .
   ```

2. Edit `compile_graspa.job` and update:
   - `--account=YOUR_ACCOUNT_HERE`: Your Quest account/project ID
   - `GRASPA_SRC_DIR`: Path to your gRASPA source directory

3. Submit the job:
   ```bash
   sbatch compile_graspa.job
   ```

4. Monitor and check results:
   ```bash
   squeue -u $USER
   cat compile_output.JOBID
   ```

#### Option B: Compile on GPU Node Directly

1. Request an interactive GPU session:
   ```bash
   srun --account=YOUR_ACCOUNT --partition=gengpu --gres=gpu:1 --time=00:30:00 --pty bash
   ```

2. Navigate to source directory and compile:
   ```bash
   cd /path/to/gRASPA/src_clean
   ./NVC_COMPILE_QUEST_VANILLA
   ```

## Verification

After compilation, verify the executable was created:
```bash
ls -lh nvc_main.x
file nvc_main.x
```

The executable should be approximately 5–6 MB in size.

## Running gRASPA

Submit a GPU job to run simulations. The binary works on **any** GPU type in the `gengpu` partition:

```bash
#!/bin/bash
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=gengpu
#SBATCH --gres=gpu:1
#SBATCH --job-name=gRASPA_run
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --error=run_error.%J
#SBATCH --output=run_output.%J

# Save SLURM GPU assignment before module purge
SAVED_CVD="${CUDA_VISIBLE_DEVICES:-}"

# Load runtime modules
module purge
module use /hpc/software/spack_v20d1/spack/share/spack/modules/linux-rhel7-x86_64/
module load nvhpc/23.3-gcc-10.4.0

# Restore GPU visibility
export CUDA_VISIBLE_DEVICES="${SAVED_CVD}"

cd /path/to/your/simulation/directory
./nvc_main.x
```

**Important:** The run script must save `CUDA_VISIBLE_DEVICES` before `module purge` and restore it afterward. Otherwise, the CUDA runtime cannot find the GPU.

## Compilation Details

### Compiler and Modules

- **NVIDIA HPC SDK**: nvhpc/23.3-gcc-10.4.0 (from spack_v20d1)
- **CUDA Toolkit**: 11.8 (bundled inside nvhpc 23.3)
- **GPU Targets**: cc80 (A100) + cc90 (H100)

### Compiler Flags

| Flag | Purpose |
|------|---------|
| `-O3` | Maximum optimization |
| `-std=c++20` | C++20 standard (for `std::filesystem` etc.) |
| `-target=gpu` | GPU code generation |
| `-gpu=cuda11.8,cc80,cc90` | CUDA 11.8 toolkit, A100 + H100 native code |
| `-Minline` | Inline function expansion |
| `-fopenmp` | OpenMP support |
| `-cuda` | CUDA support |
| `-stdpar=multicore` | Standard parallelism for multicore |

### Source Files Compiled

- `axpy.cu` — CUDA kernels for vector operations
- `main.cpp` — Main program entry point
- `read_data.cpp` — Input file reading
- `data_struct.cpp` — Data structures
- `VDW_Coulomb.cu` — Van der Waals and Coulomb interactions

### Linker Workaround

The script prepends gcc-10.3.0's `lib64/` to the linker search path. This is needed because gcc-10.4.0's `libatomic.so` (used by nvhpc/23.3) has compressed debug sections that the system linker (`ld 2.30`, RHEL 8) cannot read. The gcc-10.3.0 copy has uncompressed debug sections and works fine.

## Tested Configuration

| Component | Value |
|-----------|-------|
| **Cluster** | Quest, Northwestern University |
| **OS** | RHEL 8.10, kernel 4.18.0-553 |
| **NVIDIA Driver** | 570.86.15 (CUDA 12.8 capable) |
| **Compiler** | nvhpc/23.3-gcc-10.4.0 |
| **CUDA Toolkit** | 11.8 (bundled in nvhpc 23.3) |
| **A100 test node** | qgpu2002 — exit code 0 |
| **H100 test node** | qgpu3006 — exit code 0 |
| **Test case** | CO2 in Mg-MOF-74 (UFF), 10k cycles |
| **Test date** | April 3, 2026 |

## Troubleshooting

### "no CUDA-capable device is detected"

This means `CUDA_VISIBLE_DEVICES` or `LD_LIBRARY_PATH` was cleared by `module purge`. Fix:
```bash
# BEFORE module purge, save:
SAVED_CVD="${CUDA_VISIBLE_DEVICES:-}"

# AFTER module load, restore:
export CUDA_VISIBLE_DEVICES="${SAVED_CVD}"
```

### Linker error: "libatomic.so: file not recognized"

The compilation script already handles this. If you see this error, ensure you are using the `NVC_COMPILE_QUEST_VANILLA` script from this directory (not an older version).

### Compilation fails with Thrust error

Do NOT use nvhpc/22.7-gcc — it has an incompatible Thrust version. Use nvhpc/23.3-gcc-10.4.0.

### "no kernel image is available for execution on the device"

This means the binary was compiled without the correct compute capability for the GPU. Ensure the compile flag includes `-gpu=cc80,cc90` (both A100 and H100).

### Module not found

If `nvhpc/23.3-gcc-10.4.0` is not found:
```bash
module use /hpc/software/spack_v20d1/spack/share/spack/modules/linux-rhel7-x86_64/
module avail nvhpc
```

## Changes from Previous Quest Installation

| Change | Old (Quest/) | New (Quest_April_2026/) |
|--------|-------------|------------------------|
| Compiler | nvhpc/21.9-gcc | nvhpc/23.3-gcc-10.4.0 |
| CUDA toolkit | 11.4 (default) | 11.8 (explicit) |
| GPU targets | cc80 only (implicit) | cc80 + cc90 (explicit) |
| H100 support | No | Yes |
| libatomic fix | Not needed | gcc-10.3.0 lib64 prepended |
| CUDA_VISIBLE_DEVICES | Not handled | Save/restore around module purge |

## Additional Resources

- Quest IT Support: https://www.it.northwestern.edu/departments/it-services-support/research/computing/quest/
- gRASPA Documentation: https://zhaoli2042.github.io/gRASPA-mkdoc
- Debugging report: See `QUEST_BUILD_DEBUG.md` in the repository root

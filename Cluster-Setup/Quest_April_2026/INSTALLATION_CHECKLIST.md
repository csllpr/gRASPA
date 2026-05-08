# Installation Checklist

Use this checklist to ensure a successful gRASPA installation on Quest (April 2026).

## Pre-Installation

- [ ] Have Quest account and project ID ready
- [ ] Have access to GPU nodes on Quest (`gengpu` partition)
- [ ] Located gRASPA source code directory (`src_clean/`)

## Installation Steps

- [ ] Copied `NVC_COMPILE_QUEST_VANILLA` to `src_clean/` directory
- [ ] Made compilation script executable: `chmod +x NVC_COMPILE_QUEST_VANILLA`
- [ ] Copied `compile_graspa.job` to `src_clean/` directory
- [ ] Edited `compile_graspa.job`:
  - [ ] Updated `--account=YOUR_ACCOUNT_HERE` with your Quest account
  - [ ] Updated `GRASPA_SRC_DIR` with correct path
- [ ] Submitted compilation job: `sbatch compile_graspa.job`
- [ ] Monitored job: `squeue -u $USER`
- [ ] Checked output: `cat compile_output.JOBID`

## Verification

- [ ] Job completed successfully (check `compile_output.JOBID` for "OK nvc_main.x")
- [ ] Executable `nvc_main.x` exists in `src_clean/` directory
- [ ] Executable size is approximately 5–6 MB
- [ ] Verified executable: `file nvc_main.x` shows ELF 64-bit executable

## Post-Installation

- [ ] Created run script for gRASPA simulations
- [ ] Run script saves `CUDA_VISIBLE_DEVICES` before `module purge`
- [ ] Run script restores `CUDA_VISIBLE_DEVICES` after `module load`
- [ ] Run script loads `nvhpc/23.3-gcc-10.4.0` before executing `nvc_main.x`
- [ ] Tested running gRASPA with a sample input file on A100 or H100 node

## Troubleshooting

If compilation fails:
- [ ] Checked error log: `cat compile_error.JOBID`
- [ ] Verified module is available: `module use /hpc/software/spack_v20d1/spack/share/spack/modules/linux-rhel7-x86_64/ && module avail nvhpc`
- [ ] Confirmed job ran on a GPU node (`nvidia-smi -L` in output)
- [ ] Verified all source files are present in `src_clean/`

If runtime fails with "no CUDA-capable device":
- [ ] Confirmed `CUDA_VISIBLE_DEVICES` is saved before and restored after `module purge`
- [ ] Confirmed `nvhpc/23.3-gcc-10.4.0` is loaded in the run script
- [ ] Confirmed `--gres=gpu:1` is in the SLURM job header

## Notes

Add any custom notes or configurations here:
_____________________________________________
_____________________________________________
_____________________________________________

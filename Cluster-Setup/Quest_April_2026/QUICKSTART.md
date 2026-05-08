# Quick Start Guide: gRASPA on Quest (April 2026)

## Minimal Installation Steps

### 1. Copy Files
```bash
cd /path/to/gRASPA/src_clean
cp Cluster-Setup/Quest_April_2026/NVC_COMPILE_QUEST_VANILLA .
cp Cluster-Setup/Quest_April_2026/compile_graspa.job .
chmod +x NVC_COMPILE_QUEST_VANILLA
```

### 2. Edit Job Script
Edit `compile_graspa.job`:
- Change `YOUR_ACCOUNT_HERE` to your Quest account
- Update `GRASPA_SRC_DIR` to your source directory path

### 3. Submit Job
```bash
sbatch compile_graspa.job
```

### 4. Check Results
```bash
squeue -u $USER                    # check job status
cat compile_output.JOBID           # view output (replace JOBID)
ls -lh nvc_main.x                  # verify executable (~5-6 MB)
```

### 5. Run a Simulation
```bash
# In your run script, always do this before running nvc_main.x:
SAVED_CVD="${CUDA_VISIBLE_DEVICES:-}"
module purge
module use /hpc/software/spack_v20d1/spack/share/spack/modules/linux-rhel7-x86_64/
module load nvhpc/23.3-gcc-10.4.0
export CUDA_VISIBLE_DEVICES="${SAVED_CVD}"

./nvc_main.x
```

## That's It!

The executable `nvc_main.x` works on both **A100** and **H100** nodes. Use `--gres=gpu:1` in your SLURM scripts (no need to specify GPU type).

For detailed instructions, see [README.md](README.md).

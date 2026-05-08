# Files in Quest_April_2026 Directory

This directory contains all files needed to compile gRASPA (vanilla version) on Quest cluster, with support for both A100 and H100 GPU nodes.

## Core Files

### `NVC_COMPILE_QUEST_VANILLA`
- **Type**: Bash script
- **Purpose**: Main compilation script for gRASPA vanilla version
- **Usage**: Run directly on a GPU node or via job submission
- **Dependencies**: Requires nvhpc/23.3-gcc-10.4.0 module (spack_v20d1)
- **GPU support**: A100 (cc80) + H100 (cc90)

### `compile_graspa.job`
- **Type**: SLURM job script
- **Purpose**: Template for submitting compilation as a GPU job
- **Usage**: Edit account and paths, then `sbatch compile_graspa.job`
- **Note**: Update `YOUR_ACCOUNT_HERE` and `GRASPA_SRC_DIR` before use

## Documentation

### `README.md`
- **Type**: Markdown documentation
- **Purpose**: Comprehensive installation and usage guide
- **Contents**:
  - Prerequisites and quick start
  - Compilation details and flags
  - Runtime setup (CUDA_VISIBLE_DEVICES handling)
  - Tested configuration and troubleshooting
  - Changes from previous Quest installation

### `QUICKSTART.md`
- **Type**: Markdown documentation
- **Purpose**: Minimal quick start guide for experienced users

### `FILES.md`
- **Type**: Markdown documentation (this file)
- **Purpose**: Overview of all files in this directory

### `INSTALLATION_CHECKLIST.md`
- **Type**: Markdown documentation
- **Purpose**: Step-by-step checklist for installation verification

### `SUMMARY.md`
- **Type**: Markdown documentation
- **Purpose**: Package overview and tested configuration summary

## File Structure

```
Quest_April_2026/
├── NVC_COMPILE_QUEST_VANILLA    # Compilation script (A100 + H100)
├── compile_graspa.job            # Job submission template
├── README.md                      # Full documentation
├── QUICKSTART.md                  # Quick start guide
├── FILES.md                       # This file
├── INSTALLATION_CHECKLIST.md      # Installation checklist
└── SUMMARY.md                     # Package summary
```

## Usage Workflow

1. Copy `NVC_COMPILE_QUEST_VANILLA` to your `src_clean/` directory
2. Copy and edit `compile_graspa.job` (update account and paths)
3. Submit job: `sbatch compile_graspa.job`
4. Check results and use the generated `nvc_main.x` executable
5. In run scripts: save/restore `CUDA_VISIBLE_DEVICES` around `module purge`

For detailed instructions, see [README.md](README.md).

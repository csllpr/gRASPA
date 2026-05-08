# gRASPA Project

## Overview
gRASPA (GPU-accelerated RASPA) is a GPU-accelerated Monte Carlo simulation software for molecular adsorption in nanoporous materials (zeolites, MOFs). It is the GPU port of RASPA2 (CPU-based), built with NVIDIA's nvc++ compiler and CUDA.

## Build System
- **Compiler**: nvc++ from NVIDIA HPC SDK 24.5 at `/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/compilers/bin/nvc++`
- **Build script**: `src_clean/NVC_COMPILE` (bash script, compiles .cu and .cpp files in parallel, links to `nvc_main.x`)
- **Flags**: `-O3 -std=c++20 -target=gpu -Minline -fopenmp -cuda -stdpar=multicore`
- **Binary output**: `src_clean/nvc_main.x`
- **CUDA libs**: from HPC SDK 22.5 at `/opt/nvidia/hpc_sdk/Linux_x86_64/22.5/cuda/lib64`

## Hardware
- GPU 0: NVIDIA GeForce RTX 3080 Ti (12 GB)
- GPU 1: NVIDIA GeForce RTX 3090 (24 GB) - often loaded with batch jobs
- Use `CUDA_VISIBLE_DEVICES=0` to target the free 3080 Ti

## Source Code Structure (`src_clean/`, ~18,771 lines)
- `main.cpp` / `main.h` - Entry point, Initialize(), RunSimulation(), EndOfSimulationWrapUp()
- `read_data.cpp` / `read_data.h` - Input parsing (simulation.input, CIF, force field files)
- `data_struct.cpp` / `data_struct.h` - Core data structures (Variables, enums for MoveTypes, SIMULATION_MODE, etc.)
- `axpy.cu` / `axpy.h` - GPU kernel utilities
- `VDW_Coulomb.cu` / `VDW_Coulomb.cuh` - Van der Waals and Coulomb energy/force GPU kernels
- `Ewald_Energy_Functions.h` - Ewald summation for long-range electrostatics
- `mc_single_particle.h` - Single-particle MC moves (translation, rotation)
- `mc_swap_moves.h` / `mc_swap_utilities.h` - Insertion/deletion/swap MC moves
- `mc_cbcfc.h` - Configurational Bias / Continuous Fractional Component MC
- `mc_widom.h` - Widom test particle insertion
- `mc_box.h` - Box-related operations (Gibbs ensemble volume changes)
- `mc_utilities.h` - General MC utilities
- `lambda.h` - Lambda (alchemical) scaling for CB/CFC methods
- `fxn_main.h` - Main simulation loop functions
- `print_statistics.cuh` - Output statistics (GPU-aware)
- `write_data.h` - Output writers (restart files, LAMMPS data files)
- `equations_of_state.h` - Peng-Robinson EOS for fugacity calculations
- `TailCorrection_Energy_Functions.h` - Analytical tail corrections
- `move_struct.h` - MC move energy structures
- `maths.cuh` - Math utilities (GPU)
- `torch_allegro.h`, `cppflow_LCLin.h`, `DNN_HostGuest_Energy_Functions.h` - ML potential integration

## Key Data Structures
- `Variables` - Main simulation state container (returned by Initialize())
- `MoveTypes` enum: TRANSLATION, ROTATION, SINGLE_INSERTION, SINGLE_DELETION, INSERTION, DELETION, REINSERTION, CBCF_*, IDENTITY_SWAP, WIDOM
- `SIMULATION_MODE` enum: CREATE_MOLECULE, INITIALIZATION, EQUILIBRATION, PRODUCTION

## Simulation Input Format
gRASPA reads `simulation.input` files compatible with RASPA2 format:
```
SimulationType MonteCarlo
NumberOfCycles / NumberOfInitializationCycles / NumberOfProductionCycles
Forcefield [local | name]
UseChargesFromCIFFile yes/no
ChargeMethod Ewald
CutOff / CutOffVDW / EwaldPrecision
Framework 0: FrameworkName, UnitCells, ExternalTemperature, ExternalPressure
Component 0: MoleculeName, TranslationProbability, RotationProbability, SwapProbability, etc.
```

## Force Field Files
- `force_field.def` - Explicit pair interactions
- `force_field_mixing_rules.def` - LJ parameters per atom type + mixing rule
- `pseudo_atoms.def` - Atom type definitions (mass, charge, etc.)
- Molecule `.def` files - Molecular geometry and topology

## VDW Potential Types
### Standard 12-6 LJ (default)
```
U(r) = 4ε[(σ/r)^12 - (σ/r)^6]
```
Keyword: `lennard-jones` — used in `force_field_mixing_rules.def` and `force_field.def`

### 12-6-4 LJ (for MOFs with open metal sites)
```
U(r) = 4ε[(σ/r)^12 - (σ/r)^6] + C_4/r^4
```
Keyword: `lennard-jones-1264` — adds a charge-induced dipole r^-4 polarization term.
Based on Du, Rodriguez, Lin & Chen (J. Chem. Inf. Model. 2026).

**Implementation details:**
- C_4 is stored in the `ForceField.z` array (FFarg[2] in `maths.cuh:VDW()`)
- No GPU branching: C_4=0 for standard 12-6, non-zero for 12-6-4 (zero-cost when unused)
- Shift and tail corrections include the r^-4 contribution
- Units: ε in K, σ in Å, C_4 in K·Å⁴ (all converted internally by /1.20272430057)
- Lambda/soft-core scaling: C_4 term scales linearly with lambda (no soft-core on r^-4)

**Input format in `force_field_mixing_rules.def`:**
```
Mg  lennard-jones-1264  55.86  2.69  503.228    # epsilon sigma C_4
O   lennard-jones       48.19  3.03             # standard 12-6
```

**Pair override in `force_field.def`:**
```
Mg  O_co2  lennard-jones-1264  120.5  2.95  1250.0   # 6 tokens
Ow  Hw     lennard-jones        0.0   0.0             # 5 tokens (standard)
```

**C_4 mixing rule:** arithmetic mean: C_4,ij = (C_4,i + C_4,j) / 2

## Testing
- Examples in `Examples/` directory
- Run: `cd Examples && python run_designated_folders.py`
- Verify: `cd Examples && pytest -s` (checks energy drift < 1e-3)
- Key test cases: CO2-MFI, Methane-TMMC, Bae-Mixture, NU2000-pX-LinkerRotations, Tail-Correction, Reference_NIST_SPCE

## Related: RASPA2 (`/home/xiaoyi/RASPA2/`)
- CPU-based predecessor, version 2.0.41
- Binary: `/home/xiaoyi/RASPA2/bin/simulate`
- Contains extensive forcefield library, molecule definitions, and framework structures
- Supports MC and MD simulations
- Input format is compatible with gRASPA (gRASPA reads RASPA2-style inputs)

## Context: P_2000Pa(1).zip
A RASPA2 GCMC simulation of TIP4P/Ew water in MOF qmof-aa922ec (MgHC4O3, Mg-based MOF):
- Temperature: 298 K, Pressure: 2000 Pa
- UnitCells: 4 2 2, CutOff: 12.8 Å
- Ewald summation with precision 1e-6
- Result: ~29 mol/uc water loading (~39.7 mol/kg, ~715 mg/g)
- Uses local force field with Dreiding/UFF parameters

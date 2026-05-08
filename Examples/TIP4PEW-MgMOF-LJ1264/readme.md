# TIP4PEW Water in Mg-MOF-74 with 12-6-4 LJ Potential

## Overview

GCMC simulation of TIP4PEW water adsorption in Mg-MOF-74 (qmof-aa922ec) at 298 K / 2000 Pa,
using the **12-6-4 Lennard-Jones potential** for metal-adsorbate interactions at open metal sites.

The 12-6-4 potential extends the standard 12-6 LJ by adding a charge-induced dipole term (r^-4),
capturing polarization effects between the charged open metal site and guest molecules:

```
E(r) = 4epsilon * [(sigma/r)^12 - (sigma/r)^6] - C_4/r^4 + q_i*q_j/(4*pi*eps_0*r)
```

This example validates gRASPA's 12-6-4 implementation against RASPA2 v2.0.47 and experimental
water adsorption data in Mg-MOF-74 (see Figure 7b in Du et al.).

## Reference

**Du, M.; Rodriguez, A.; Lin, M. Z.; Chen, H.**
*A Transferable Force Field for Simulating Adsorption in Metal-Organic Frameworks with Open Metal Sites Based on the 12-6-4 Lennard-Jones Potential.*
J. Chem. Inf. Model. **2026**, 66, 1704-1714.
DOI: [10.1021/acs.jcim.5c02893](https://doi.org/10.1021/acs.jcim.5c02893)

### Key findings from the paper:
- The 12-6-4 FF is parametrized against DFT-derived potential energy surfaces (PES) for
  60 metal-guest combinations (12 metals x 5 guests) using particle swarm optimization (PSO)
- Compared to DREIDING/UFF, the 12-6-4 FF reduces mean absolute errors (MAEs) from >100 kJ/mol
  to <10 kJ/mol for OMS-guest binding energies
- The C_4 parameter in GCMC simulations is charge-scaled:
  `C_4^MOF = C_4^TMC * (q_M^MOF / q_M^TMC)^2` to account for different partial charges
  between the training cluster (TMC) and the actual MOF
- Water adsorption isotherms in Mg-MOF-74 (Figure 7b) show excellent agreement with
  experimental data from Glover et al. and Schoenecker et al.

### RASPA2 implementation:
- Modified RASPA2 source: https://github.com/haoyuanchen/RASPA-tools/tree/master/LJ1264Potential

## 12-6-4 FF Parameters (Table 1 in Du et al.)

Optimized parameters for Mg(II)-guest atom pairs from the paper:

| Metal | Guest atom | C_12 [K A^12] | C_6 [K A^6] | C_4 [K A^4] |
|-------|-----------|---------------|-------------|-------------|
| Mg(II) | O (guest) | 10311.42 | 47.51 | 206.60 |
| Mg(II) | C (guest) | 18786.72 | 86.84 | 319.76 |

These are the TMC-fitted values. For GCMC in a specific MOF, C_4 is charge-scaled by
(q_Mg^MOF / q_Mg^TMC)^2, where q_Mg^TMC = +2 (formal oxidation state).

## Input Setup

### simulation.input
- `UseLJ1264 yes` enables the polynomial VDW mode in gRASPA
- Temperature: 298 K, Pressure: 2000 Pa
- Framework: qmof-aa922ec (Mg-MOF-74), UnitCells: 4 2 2
- CutOff VDW: 12.8 A, CutOff Coulomb: 12.0 A, Ewald precision: 1e-6
- 10,000 initialization + 100,000 production cycles

### force_field.def (GENERIC2_HC pair override)
The Mg-Ow interaction is overridden with the 12-6-4 parameters in RASPA2's GENERIC2_HC format:

```
Mg  OwH2O_TIP4PEW  GENERIC2_HC  0  0  5458  2234  0  -14427061
```

GENERIC2_HC formula: `U(r) = p0*exp(-p1*r) - p2/r^4 - p3/r^6 - p4/r^8 - p5/r^12`
- p2 = 5458 [K A^4] — C_4 charge-induced dipole term
- p3 = 2234 [K A^6] — C_6 dispersion term
- p5 = -14427061 [K A^12] — C_12 repulsion (negative sign convention in GENERIC2_HC)
- p0, p1, p4 = 0 (unused exponential and r^-8 terms)

### force_field_mixing_rules.def
Standard LJ parameters for all atom types (UFF for Mg, Dreiding for framework H/C/O,
TIP4PEW for water). Truncated at cutoff, no tail corrections.

### pseudo_atoms.def
Framework charges from DDEC analysis of the CIF file (`UseChargesFromCIFFile yes`).
TIP4PEW water: Hw = +0.52422e, Mw = -1.04844e (M-site model).

## Results: gRASPA vs RASPA2 vs Paper

All simulations at 298 K, 2000 Pa, CutOff VDW 12.8 A, truncated, no tail corrections.

### Loading comparison

| Property | gRASPA (100k prod) | gRASPA (20k prod) | RASPA2 v2.0.47 (100k) | Paper Fig 7b |
|----------|--------------------|--------------------|------------------------|--------------|
| Loading (molecules) | 473.9 +/- 13.9 | 460.1 +/- 15.1 | 487.4 +/- 10.4 | - |
| Loading (mol/kg) | 40.7 +/- 1.2 | 39.5 +/- 1.3 | 41.8 +/- 0.9 | ~35-42 |
| Loading (mg/g) | 732.3 +/- 21.5 | 711.0 +/- 23.3 | 753.1 +/- 16.1 | - |

All runs are mutually consistent within error bars. The paper (Figure 7b) shows H2O
adsorption in Mg-MOF-74 at 2000 Pa with the 12-6-4 FF gives ~35-42 mol/kg, matching
both gRASPA and RASPA2 results. Experimental data from Glover et al. and Schoenecker
et al. also fall in this range at saturation conditions.

### Energy comparison (100k production runs)

| Property | gRASPA | RASPA2 v2.0.47 | Diff |
|----------|--------|----------------|------|
| GG VDW (K) | +524,906 | +521,432 | 0.7% |
| GG Coulomb (K) | -2,924,905 | -2,932,894 | 0.3% |
| Energy drift (K) | 0.0 | ~5.6e-9 | Both excellent |

### Validation notes
- 12-6-4 formula verified to machine precision (<10^-13 K) against haoyuanchen potentials.c
- GG VDW close match confirms standard LJ mixing rules are correct
- Loading overlap within error bars confirms physical consistency across gRASPA, RASPA2, and the paper
- HG energies differ slightly due to different Ewald implementations and stochastic sampling
- The 20k production test run (1451 s on RTX 3080 Ti) reproduces the 100k result within error bars

## Files

| File | Description |
|------|-------------|
| `simulation.input` | Main configuration (UseLJ1264 yes) |
| `force_field.def` | GENERIC2_HC pair override for Mg-Ow |
| `force_field_mixing_rules.def` | Standard LJ parameters (truncated, no tail corrections) |
| `pseudo_atoms.def` | Atom types (charges from CIF via DDEC) |
| `qmof-aa922ec.cif` | Framework structure with DDEC charges |
| `TIP4PEW.def` | TIP4P/Ew water molecule definition |
| `output.txt` | gRASPA output (100k production cycles) |
| `RASPA2-reference-output.txt` | RASPA2 v2.0.47 output (same conditions) |

## How to Run

```bash
cd Examples/TIP4PEW-MgMOF-LJ1264
CUDA_VISIBLE_DEVICES=0 ../../src_clean/nvc_main.x
```

Expected runtime: ~90 minutes on RTX 3080 Ti (100k production cycles).

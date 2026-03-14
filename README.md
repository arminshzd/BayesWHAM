# BayesWHAM

A Bayesian implementation of the Weighted Histogram Analysis Method (WHAM) for estimating multidimensional molecular free energy surfaces from umbrella sampling simulations.

> **This repository is a Python 3 refactor of the original BayesWHAM code written by Andrew L. Ferguson (University of Illinois at Urbana-Champaign, 2016).** The original Python 2.7 scripts are preserved in `src/python2/`. The refactored package in `bayeswham/` is algorithmically identical — same WHAM iteration, Metropolis-Hastings sampler, prior implementations, and output formats — but is restructured as an installable Python package with YAML-based configuration and improved diagnostics. All scientific credit belongs to the original author.

---

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [BayesWHAM Usage](#bayeswham-usage)
  - [YAML Configuration (Recommended)](#yaml-configuration-recommended)
  - [Legacy CLI](#legacy-cli)
  - [Output Files](#bayeswham-output-files)
- [BayesReweight Usage](#bayesreweight-usage)
  - [YAML Configuration](#yaml-configuration-1)
  - [Legacy CLI](#legacy-cli-1)
  - [Output Files](#bayesreweight-output-files)
- [Plotting](#plotting)
- [PLUMED Integration](#plumed-integration)
- [Input File Formats](#input-file-formats)
- [Energy Units](#energy-units)
- [Examples](#examples)
- [Citations](#citations)
- [License](#license)

---

## Overview

**BayesWHAM** provides statistically optimal estimates of free energy surfaces (FES) from molecular dynamics umbrella sampling data. It solves the generalized WHAM self-consistency equations under a Bayesian prior and estimates uncertainties by Metropolis-Hastings (MH) sampling from the posterior distribution.

The package provides five command-line tools:

| Tool | Description |
|---|---|
| `bayeswham` | Core WHAM analysis: computes MAP estimate and MH uncertainty samples of the FES |
| `bayesreweight` | Projects the FES onto arbitrary variables beyond the umbrella sampling coordinates |
| `bayeswham-plot` | Plotting utility for 1D, 2D, and 3D BayesWHAM results |
| `bayesreweight-plot` | Plotting utility for reweighted projection results |
| `bayeswham-plumed` | Converts PLUMED COLVAR files directly into bayeswham input files |

---

## Features

- n-dimensional free energy estimation
- Native support for periodic and non-periodic coordinates
- Bayesian priors: uniform (`none`), symmetric Dirichlet, and Gaussian
- Rigorous uncertainty quantification via Metropolis-Hastings posterior sampling
- YAML-based configuration (recommended) with legacy CLI backward compatibility
- Direct PLUMED COLVAR ingestion via `bayeswham-plumed`
- Python 3.7+ compatible, installable as a package

---

## Installation

### From source

```bash
git clone https://github.com/arminshzd/BayesWHAM.git
cd BayesWHAM

# Core install
pip install -e .

# With plotting support (matplotlib)
pip install -e ".[plotting]"

# All optional dependencies (plotting + dev/test)
pip install -e ".[all]"
```

### Requirements

| Package | Version | Notes |
|---|---|---|
| Python | >= 3.7 | |
| NumPy | >= 1.19.0 | |
| SciPy | >= 1.5.0 | Required for connectivity analysis |
| PyYAML | >= 5.3.0 | |
| Matplotlib | >= 3.3.0 | Optional, for plotting tools |

---

## Quick Start

```bash
# 1. Run BayesWHAM to compute the free energy surface
bayeswham --config config.yaml

# 2. (Optional) Reweight onto different coordinates
bayesreweight --config config_reweight.yaml

# 3. Plot results
bayeswham-plot
bayesreweight-plot
```

---

## BayesWHAM Usage

### YAML Configuration (Recommended)

Create a `config.yaml` file:

```yaml
umbrella_sampling:
  dim: 1                        # Dimensionality of umbrella sampling
  periodicity: [1]              # Per-dimension: 1 = periodic, 0 = non-periodic
  temperature: 298.0            # Temperature in Kelvin
  # boltzmann_constant: 0.0083144621  # Optional; default is kJ/mol·K

input_files:
  harmonic_biases: "./data/harmonic_biases.txt"
  hist_bin_edges: "./data/hist_binEdges.txt"
  hist_dir: "./data/hist"

wham:
  tolerance: 1.0e-15            # Convergence tolerance for WHAM self-consistency
  max_iterations: 1000000       # Maximum WHAM iterations

metropolis_hastings:
  steps: 1000000                # Total MH sampling steps
  save_modulus: 1000            # Save every Nth sample
  print_modulus: 1000           # Print progress every Nth step
  max_step_size: 5.0e-4         # Max step size in probability space
  random_seed: 200184           # RNG seed for reproducibility

prior:
  type: "Dirichlet"             # 'none' (uniform), 'Dirichlet', or 'Gaussian'
  alpha: 2.0                    # Hyperparameter: concentration (Dirichlet) or scale (Gaussian)

# output_dir: "./output"        # Optional; defaults to current directory
```

Run:

```bash
bayeswham --config config.yaml
```

### Legacy CLI

```bash
bayeswham <dim> <periodicity> <T> <harmonicBiasesFile> <histBinEdgesFile> <histDir> \
          <tol_WHAM> <maxIter_WHAM> <steps_MH> <saveMod_MH> <printMod_MH> \
          <maxDpStep_MH> <seed_MH> <prior> <alpha> [kB]
```

**Arguments:**

| # | Name | Description |
|---|---|---|
| 1 | `dim` | Integer dimensionality of umbrella sampling |
| 2 | `periodicity` | List, e.g. `[1]` or `[0,1]` — 1 = periodic, 0 = non-periodic per dimension |
| 3 | `T` | Temperature in Kelvin |
| 4 | `harmonicBiasesFile` | Path to harmonic bias parameters file |
| 5 | `histBinEdgesFile` | Path to histogram bin edges file |
| 6 | `histDir` | Directory containing histogram files (`hist_1.txt`, `hist_2.txt`, ...) |
| 7 | `tol_WHAM` | WHAM convergence tolerance (e.g., `1e-15`) |
| 8 | `maxIter_WHAM` | Maximum WHAM iterations (e.g., `1000000`) |
| 9 | `steps_MH` | Number of MH sampling steps |
| 10 | `saveMod_MH` | Save every Nth MH sample |
| 11 | `printMod_MH` | Print progress every Nth step |
| 12 | `maxDpStep_MH` | Maximum step size in probability space (e.g., `5e-4`) |
| 13 | `seed_MH` | Integer random seed |
| 14 | `prior` | `'none'`, `'Dirichlet'`, or `'Gaussian'` |
| 15 | `alpha` | Prior hyperparameter |
| 16 | `kB` | *(Optional)* Boltzmann constant; default `0.0083144621` kJ/mol·K |

### BayesWHAM Output Files

All files are written to the current directory (or `output_dir` if specified):

| File | Description |
|---|---|
| `hist_binCenters.txt` | Bin center coordinates for each dimension |
| `hist_binWidths.txt` | Bin widths for each dimension |
| `p_MAP.txt` | MAP estimate of the probability distribution |
| `pdf_MAP.txt` | MAP estimate of the probability density function |
| `betaF_MAP.txt` | MAP estimate of the reduced free energy (βF = F/kT) |
| `f_MAP.txt` | MAP estimates of partition function ratios |
| `logL_MAP.txt` | Log-likelihood at the MAP estimate |
| `p_MH.txt` | MH samples of probability distribution |
| `pdf_MH.txt` | MH samples of probability density function |
| `betaF_MH.txt` | MH samples of reduced free energy |
| `f_MH.txt` | MH samples of partition function ratios |
| `logL_MH.txt` | Log-likelihood trajectory over MH samples |
| `step_MH.txt` | MH step indices for saved samples |

---

## BayesReweight Usage

`bayesreweight` takes the partition function ratios from `bayeswham` output and reweights the biased simulation data onto an arbitrary set of projection coordinates.

### YAML Configuration

Create a `config_reweight.yaml` file:

```yaml
temperature: 298.0
# boltzmann_constant: 0.0083144621  # Optional

umbrella:
  dim: 1
  periodicity: [1]
  files:
    harmonic_biases: "./data/harmonic_biases.txt"
    trajectory_dir: "./data/traj"          # Trajectory files: traj_1.txt, traj_2.txt, ...
    hist_bin_edges: "./data/hist_binEdges.txt"
    hist_dir: "./data/hist"
    f_map: "./output/f_MAP.txt"            # From bayeswham output
    f_mh: "./output/f_MH.txt"             # From bayeswham output

projection:
  files:
    trajectory_dir: "./data/traj_projection"   # Projection variable trajectories
    hist_bin_edges: "./data/hist_binEdges_projection.txt"

# output_dir: "./output"
```

Run:

```bash
bayesreweight --config config_reweight.yaml
```

### Legacy CLI

```bash
bayesreweight <T> <dim_UMB> <periodicity_UMB> <harmonicBiasesFile_UMB> <trajDir_UMB> \
              <histBinEdgesFile_UMB> <histDir_UMB> <fMAPFile_UMB> <fMHFile_UMB> \
              <trajDir_PROJ> <histBinEdgesFile_PROJ>
```

**Arguments:**

| # | Name | Description |
|---|---|---|
| 1 | `T` | Temperature in Kelvin |
| 2 | `dim_UMB` | Dimensionality of umbrella sampling coordinates |
| 3 | `periodicity_UMB` | Periodicity per umbrella dimension |
| 4 | `harmonicBiasesFile_UMB` | Harmonic bias parameters (same format as BayesWHAM) |
| 5 | `trajDir_UMB` | Directory with umbrella trajectory files (`traj_1.txt`, ...) |
| 6 | `histBinEdgesFile_UMB` | Umbrella histogram bin edges |
| 7 | `histDir_UMB` | Directory with umbrella histograms |
| 8 | `fMAPFile_UMB` | Path to `f_MAP.txt` from BayesWHAM output |
| 9 | `fMHFile_UMB` | Path to `f_MH.txt` from BayesWHAM output |
| 10 | `trajDir_PROJ` | Directory with projection variable trajectories (must match length of umbrella trajectories) |
| 11 | `histBinEdgesFile_PROJ` | Bin edges for the projection histogram |

### BayesReweight Output Files

| File | Description |
|---|---|
| `hist_binCenters_PROJ.txt` | Projection bin center coordinates |
| `hist_binWidths_PROJ.txt` | Projection bin widths |
| `p_PROJ_MAP.txt` | Reweighted MAP probability distribution |
| `pdf_PROJ_MAP.txt` | Reweighted MAP probability density |
| `betaF_PROJ_MAP.txt` | Reweighted MAP reduced free energy |
| `p_PROJ_MH.txt` | Reweighted MH probability samples |
| `pdf_PROJ_MH.txt` | Reweighted MH probability density samples |
| `betaF_PROJ_MH.txt` | Reweighted MH free energy samples |

---

## Plotting

```bash
# Plot BayesWHAM results (reads output files from current directory by default)
bayeswham-plot

# Plot BayesReweight results
bayesreweight-plot
```

Both plotters support 1D, 2D, and 3D visualizations and output `.jpg` and `.eps` files. They read from default filenames in the current directory; pass custom paths as arguments if needed.

---

## PLUMED Integration

`bayeswham-plumed` converts PLUMED COLVAR files from umbrella sampling runs directly into all input files needed by `bayeswham`, and writes a ready-to-use `bayeswham_config.yaml`.

### Input file format

Create a plain text file listing one umbrella window per line:

```
# colvar_path                         umb_center   umb_k
/path/to/window_01/COLVAR             -170.0       0.0239
/path/to/window_02/COLVAR             -150.0       0.0239
/path/to/window_03/COLVAR             -130.0       0.0239
...
```

Lines beginning with `#` are ignored. Fields are whitespace-delimited.

### Usage

```bash
bayeswham-plumed windows.txt [options]
```

**Options:**

| Flag | Default | Description |
|---|---|---|
| `--variable NAME` | `cv1` | Column name to extract from the `#! FIELDS` header |
| `--output-dir DIR` | `./bayeswham_input` | Directory to write all output files |
| `--nbins N` | `50` | Number of histogram bins |
| `--temperature T` | `300.0` | Temperature in Kelvin |
| `--skip N` | `0` | Number of initial frames to skip |
| `--discard-fraction F` | `0.1` | Fraction of the leading trajectory to discard for equilibration (applied after `--skip`, always from the start) |
| `--bin-range MIN MAX` | auto | Manual histogram range; auto-detected from data if omitted |

> **Force constant units (`umb_k`):** Must be in **energy_unit / CV_unit²**, matching the Boltzmann constant used by `bayeswham`. With the default `kB = 0.0083144621 kJ/mol·K`, `umb_k` should be in **kJ/mol per (CV unit)²** — consistent with PLUMED's `RESTRAINT` `kappa` output. If your simulation used kcal/mol, set `boltzmann_constant: 0.001987204` in the generated `bayeswham_config.yaml`.

> **Equilibration discard (`--discard-fraction`):** Removes the specified fraction of frames from the **beginning** of each trajectory (after any `--skip` frames). For example, `--discard-fraction 0.2` discards the first 20% of each window's trajectory. There is no end-trimming; all retained frames run to the end of the COLVAR file.

### Example

```bash
bayeswham-plumed windows.txt \
    --variable cv1 \
    --nbins 60 \
    --temperature 300 \
    --discard-fraction 0.2 \
    --output-dir ./wham_input

# Then run bayeswham with the generated config
bayeswham --config ./wham_input/bayeswham_config.yaml
```

### Outputs

`bayeswham-plumed` writes the following into `--output-dir`:

| Path | Description |
|---|---|
| `hist_binEdges.txt` | Histogram bin edges |
| `bias/harmonic_biases.txt` | Umbrella centers and force constants |
| `hist/hist_1.txt`, `hist_2.txt`, ... | Per-window histogram counts |
| `bayeswham_config.yaml` | Ready-to-use bayeswham config (edit prior/MH settings as needed) |
| `simulation_summary.txt` | Human-readable summary of loaded windows |

The generated `bayeswham_config.yaml` uses sensible defaults (Dirichlet prior, α=2, 500k MH steps). Review and adjust the `metropolis_hastings` and `prior` sections before running `bayeswham` on production data.

### Using bayesreweight after bayeswham-plumed

`bayeswham-plumed` produces everything needed to run `bayeswham` and recover F in the umbrella sampling coordinate. If you also want to project the FES onto **different coordinates** using `bayesreweight`, you will need single-column trajectory files for each umbrella window in the projection variable(s).

These can be extracted directly from your COLVAR files. For example, to extract column `cv2` from each COLVAR:

```bash
# Get the column index of cv2 from the header (0-indexed, excluding the #! FIELDS token)
# Example header: #! FIELDS time cv1 cv2 ...
#   time=col1, cv1=col2, cv2=col3 → awk field $3

mkdir -p traj_cv2
for i in $(seq 1 N); do
    awk 'NR>1 && !/^#/ {print $3}' /path/to/window_${i}/COLVAR > traj_cv2/traj_${i}.txt
done
```

Adjust the awk field index (`$3`) to match the column position of your projection variable in the COLVAR header. The trajectory files must be in the same order and have the same number of rows as the umbrella windows passed to `bayeswham-plumed` (after applying the same `--skip` and `--discard-fraction`).

---

## Input File Formats

### `harmonic_biases.txt`

Space-delimited, one umbrella window per line:

```
simulation_ID  center_1 [center_2 ...]  force_const_1 [force_const_2 ...]
```

Example (1D):
```
1  -170.0  0.0239
2  -150.0  0.0239
...
```

Force constants must match the energy units of your Boltzmann constant (default: kJ/mol·deg²).

### `hist_binEdges.txt`

One row per dimension, space-delimited bin edge values. For M bins in dimension k, provide M+1 edge values:

```
-180.0 -162.0 -144.0 ... 180.0
```

### Histogram files (`hist_dir/hist_1.txt`, `hist_2.txt`, ...)

Row vectors in row-major order (last index changes fastest for multi-dimensional histograms). One file per umbrella window.

### Trajectory files (`traj_dir/traj_1.txt`, `traj_2.txt`, ...)

N rows × dim columns, space-delimited coordinate values. One file per umbrella window. Projection trajectories must have the same number of rows as the corresponding umbrella trajectory files.

---

## Energy Units

The default Boltzmann constant is `kB = 0.0083144621 kJ/mol·K`. Override it to use different energy units:

| Value | Units |
|---|---|
| `0.0083144621` | kJ/mol·K (default) |
| `0.001987204` | kcal/mol·K |
| `1.380649e-23` | J/K |

In YAML:
```yaml
umbrella_sampling:
  temperature: 298.0
  boltzmann_constant: 0.001987204
```

In legacy CLI, pass `kB` as the optional 16th argument to `bayeswham`.

Ensure that force constants in `harmonic_biases.txt` use consistent units.

---


## Citations

If you use this software, please cite the original BayesWHAM paper:

> A.L. Ferguson, "BayesWHAM: A Bayesian approach for free energy estimation, reweighting, and uncertainty quantification in the weighted histogram analysis method," *J. Comput. Chem.* **38**(18), 1583–1605 (2017).

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for the full text.

```
The MIT License (MIT)

Copyright (c) 2016 Andrew L Ferguson

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
```

---

## Version History

- **v1.1.0** (2025) — Python 3 refactor by Armin Shayesteh Zadeh: installable package structure, YAML configuration, `main()` entry points, enhanced diagnostics. Original algorithms unchanged.
- **v1.0** (2016) — Initial release by Andrew L. Ferguson (Python 2.7)

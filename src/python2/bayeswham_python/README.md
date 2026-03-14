# README #



### Synopsis ###

**BayesWHAM.py** - a Bayesian implementation of the weighted histogram analysis method (WHAM) to make statistically optimal estimates of multidimensional molecular free energy surface (FES) from an ensemble of umbrella sampling simulations. The code functions in n-dimensions, natively supports uniform, symmetric Dirichlet, and Gaussian priors, and rigorously estimates uncertainties by Metropolis-Hastings sampling of the Bayes posterior. The source code may be straightforwardly modified to support arbitrary user-defined priors by augmenting two switch statements defining the MAP expression and the Metropolis acceptance criterion.

**BayesWHAM_plotter.py** - a utility code to plot the 1, 2, and 3-dimensional output of BayesWHAM.py

**BayesReweight.py** - a tool to project the unbiased FES calculated by BayesWHAM into any number of arbitrary variables beyond those in which umbrella sampling was conducted.

**BayesReweight_plotter.py** - a utility code to plot the 1, 2, and 3-dimensional output of BayesReweight.py



### Version ###

v1.0 - 08/2016 (Updated for Python 3 compatibility - 2025)



### Installation ###

No installation required. These are Python 3 scripts that can be called and executed directly from the Python interpreter.

**Requirements:**
- Python 3.x or higher
- NumPy
- SciPy (for BayesWHAM.py connectivity analysis)
- Matplotlib (for plotter scripts)

**Installation via conda/pip:**
```bash
# Using conda
conda install numpy scipy matplotlib

# Or using pip
pip install numpy scipy matplotlib
```



### I/O Requirements ###

#### BayesWHAM.py - Input Files

**Required inputs (15 command-line arguments):**

1. **dim** - Dimensionality of umbrella sampling (integer)
2. **periodicity** - Periodicity in each dimension, e.g., `[1]` for 1D periodic or `[0,1]` for 2D with second dimension periodic
3. **T** - Temperature in Kelvin (float)
4. **harmonicBiasesFile** - Path to text file with umbrella simulation parameters
   - Format: Space-delimited, one simulation per line
   - Columns: `simulation_ID  center_1 ... center_dim  force_const_1 ... force_const_dim`
   - Force constants in kJ/mol·(units)²
5. **histBinEdgesFile** - Path to histogram bin edge definitions
   - Format: One row per dimension, space-delimited bin edge values
   - For M bins in dimension k, provide (M+1) edge values
6. **histDir** - Directory containing histogram files `hist_1.txt`, `hist_2.txt`, etc.
   - Format: Row vectors in row-major order (last index changes fastest)
7. **tol_WHAM** - Convergence tolerance for WHAM equations (e.g., 1e-15)
8. **maxIter_WHAM** - Maximum iterations for WHAM convergence (e.g., 1000000)
9. **steps_MH** - Number of Metropolis-Hastings sampling steps (e.g., 1000000)
10. **saveMod_MH** - Save every Nth MH sample (e.g., 1000)
11. **printMod_MH** - Print progress every Nth MH step (e.g., 1000)
12. **maxDpStep_MH** - Maximum step size in probability space (e.g., 5e-4)
13. **seed_MH** - Random seed for MH sampling (integer)
14. **prior** - Prior type: `'none'`, `'Dirichlet'`, or `'Gaussian'`
15. **alpha** - Hyperparameter for prior distribution (float)

#### BayesWHAM.py - Output Files

All output files are written to the current working directory:

- **hist_binCenters.txt** - Bin center coordinates for each dimension
- **hist_binWidths.txt** - Bin widths for each dimension
- **p_MAP.txt** - MAP estimate of probability distribution
- **pdf_MAP.txt** - MAP estimate of probability density function
- **betaF_MAP.txt** - MAP estimate of free energy (βF)
- **f_MAP.txt** - MAP estimates of partition function ratios
- **logL_MAP.txt** - Log-likelihood at MAP estimate
- **p_MH.txt** - Metropolis-Hastings samples of probability distribution
- **pdf_MH.txt** - MH samples of probability density function
- **betaF_MH.txt** - MH samples of free energy
- **f_MH.txt** - MH samples of partition function ratios
- **logL_MH.txt** - Log-likelihood trajectory over MH samples
- **step_MH.txt** - MH step numbers for saved samples

#### BayesReweight.py - Input Files

**Required inputs (11 command-line arguments):**

1. **T** - Temperature in Kelvin
2. **dim_UMB** - Dimensionality of umbrella sampling coordinates
3. **periodicity_UMB** - Periodicity in each umbrella dimension
4. **harmonicBiasesFile_UMB** - Umbrella bias parameters (same format as BayesWHAM)
5. **trajDir_UMB** - Directory with umbrella trajectory files `traj_1.txt`, `traj_2.txt`, etc.
   - Format: N rows × dim_UMB columns, space-delimited
6. **histBinEdgesFile_UMB** - Umbrella histogram bin edges
7. **histDir_UMB** - Directory with umbrella histograms
8. **fMAPFile_UMB** - Path to `f_MAP.txt` from BayesWHAM output
9. **fMHFile_UMB** - Path to `f_MH.txt` from BayesWHAM output
10. **trajDir_PROJ** - Directory with projection variable trajectories
    - Format: N rows × dim_PROJ columns (must match trajectory length in trajDir_UMB)
11. **histBinEdgesFile_PROJ** - Bin edges for projection histogram

#### BayesReweight.py - Output Files

- **hist_binCenters_PROJ.txt** - Projection bin centers
- **hist_binWidths_PROJ.txt** - Projection bin widths
- **p_PROJ_MAP.txt** - Reweighted MAP probability distribution
- **pdf_PROJ_MAP.txt** - Reweighted MAP probability density
- **betaF_PROJ_MAP.txt** - Reweighted MAP free energy
- **p_PROJ_MH.txt** - Reweighted MH probability samples
- **pdf_PROJ_MH.txt** - Reweighted MH probability density samples
- **betaF_PROJ_MH.txt** - Reweighted MH free energy samples

#### Plotter Scripts

Both plotter scripts accept either command-line arguments or use default filenames:

**BayesWHAM_plotter.py** - Plots MAP estimates and uncertainty quantification
- Default inputs: `hist_binCenters.txt`, `pdf_MAP.txt`, `betaF_MAP.txt`, `f_MAP.txt`, etc.
- Outputs: Multiple `.jpg` and `.eps` plots for 1D, 2D, or 3D visualizations

**BayesReweight_plotter.py** - Plots reweighted projections
- Default inputs: `hist_binCenters_PROJ.txt`, `pdf_PROJ_MAP.txt`, `betaF_PROJ_MAP.txt`, etc.
- Outputs: Projection plots in `.jpg` and `.eps` formats

### Detailed I/O Specifications ###

Full details of file formats, data structures, and formatting requirements are provided in the docstring headers of each script. Run any script without arguments to see usage information.



### YAML Configuration (Recommended) ###

Both **BayesWHAM.py** and **BayesReweight.py** support YAML configuration files for easier parameter management. This is the recommended approach as it provides better organization and readability compared to long command-line arguments.

#### BayesWHAM.py Configuration

Create a YAML configuration file (e.g., `config.yaml`) with the following structure:

```yaml
# BayesWHAM Configuration File
umbrella_sampling:
  dim: 1                        # Dimensionality of umbrella sampling
  periodicity: [1]              # Periodicity in each dimension (1=periodic, 0=non-periodic)
  temperature: 298.0            # Temperature in Kelvin

input_files:
  harmonic_biases: "./diala_phi_EXAMPLE/bias/harmonic_biases.txt"
  hist_bin_edges: "./diala_phi_EXAMPLE/hist/hist_binEdges.txt"
  hist_dir: "./diala_phi_EXAMPLE/hist"

wham:
  tolerance: 1.0e-15           # Convergence tolerance for WHAM equations
  max_iterations: 1000000      # Maximum iterations for WHAM convergence

metropolis_hastings:
  steps: 1000000               # Number of MH sampling steps
  save_modulus: 1000           # Save every Nth sample
  print_modulus: 1000          # Print progress every Nth step
  max_step_size: 5.0e-4        # Maximum step size in probability space
  random_seed: 200184          # Random seed for reproducibility

prior:
  type: "Dirichlet"            # Prior type: 'none', 'Dirichlet', or 'Gaussian'
  alpha: 2.0                   # Hyperparameter (concentration for Dirichlet, scale for Gaussian)

# Optional: Specify output directory (defaults to current directory)
# output_dir: "./output"
```

**Execution:**
```bash
python3 BayesWHAM.py --config config.yaml
```

#### BayesReweight.py Configuration

Create a YAML configuration file (e.g., `config_reweight.yaml`):

```yaml
# BayesReweight Configuration File
# Temperature
temperature: 298.0             # Temperature in Kelvin

# Umbrella sampling coordinates
umbrella:
  dim: 1                       # Dimensionality of umbrella sampling
  periodicity: [1]             # Periodicity in each dimension

  # Input files for umbrella sampling
  files:
    harmonic_biases: "./diala_phi_EXAMPLE/bias/harmonic_biases.txt"
    trajectory_dir: "./diala_phi_EXAMPLE/traj"
    hist_bin_edges: "./diala_phi_EXAMPLE/hist/hist_binEdges.txt"
    hist_dir: "./diala_phi_EXAMPLE/hist"
    f_map: "./diala_phi_EXAMPLE/BayesWHAM_OUTPUT/f_MAP.txt"
    f_mh: "./diala_phi_EXAMPLE/BayesWHAM_OUTPUT/f_MH.txt"

# Projection coordinates
projection:
  # Input files for projection
  files:
    trajectory_dir: "./diala_phi_EXAMPLE/traj_aux__phi_psi"
    hist_bin_edges: "./diala_phi_EXAMPLE/hist_binEdges_phi_psi.txt"

# Optional: Specify output directory (defaults to current directory)
# output_dir: "./output"
```

**Execution:**
```bash
python3 BayesReweight.py --config config_reweight.yaml
```

#### Legacy Command-Line Interface

Both scripts maintain backward compatibility with the original command-line argument interface. See the examples below for the legacy CLI syntax.



### Example 1 - alanine dipeptide 1D umbrella sampling in phi dihedral###

* **Umbrella sampling input data**

*Ensemble of 18 1D umbrella sampling calculations in phi dihedral*

./diala_phi_EXAMPLE/bias/harmonic_biases.txt - list of umbrella simulation ID, location of harmonic biasing center (degrees), and harmonic force constant (kJ/mol.degrees^2), one entry per line

./diala_phi_EXAMPLE/traj - trajectories of umbrella sampling variable phi in each umbrella sampling calculation

./diala_phi_EXAMPLE/hist/hist_binEdges.txt - bin edges of rectilinear grid used to create phi histograms over each umbrella sampling trajectory

./diala_phi_EXAMPLE/hist - phi histograms compiled over phi trajectories in ./diala_phi_EXAMPLE/traj

./diala_phi_EXAMPLE/traj_aux__phi_psi - trajectories of projection variables [phi,psi] into which we wish to reweight the FES recorded over the course of each umbrella sampling calculation at the same frequency as those in ./diala_phi_EXAMPLE/traj; the projection variables are arbitrary and need not contain the umbrella sampling variable


* **BayesWHAM.py**

*Calculation of 1D FES F(phi) for alanine dipeptide from 18 1D umbrella sampling simulations at values of phi=[-170:20:170]. Phi is a backbone dihedral angle with range [-180 deg, 180 deg] and 360 deg periodicity.* 

*We solve the generalized WHAM equations to a tolerance of 1E-15 under a Dirichlet prior with concentration parameter alpha = 2 (i.e., adding 1 pseudo count to each bin). Uncertainties are estimated by performing 1E6 rounds of Metropolis-Hastings sampling from the Bayes posterior and saving every 1E3 realizations. In practice, many more rounds of sampling should be performed until the log likelihood stabilizes indicating the Markov chain burn-in period has passed and samples are being drawn from the stationary distribution.*

*Copies of all files generated by this calculation are provided in ./diala_phi_EXAMPLE/BayesWHAM_OUTPUT*

***Execution:***

dim=1;

periodicity=[1];

T=298;

harmonicBiasesFile='./diala_phi_EXAMPLE/bias/harmonic_biases.txt';

histBinEdgesFile='./diala_phi_EXAMPLE/hist/hist_binEdges.txt';

histDir='./diala_phi_EXAMPLE/hist';

tol_WHAM=1E-15;

maxIter_WHAM=1E6;

steps_MH=1E6;

saveMod_MH=1E3;

printMod_MH=1E3;

maxDpStep_MH=5E-4;

seed_MH=200184;

prior='Dirichlet';

alpha=2;

python3 BayesWHAM.py $dim $periodicity $T $harmonicBiasesFile $histBinEdgesFile $histDir $tol_WHAM $maxIter_WHAM $steps_MH $saveMod_MH $printMod_MH $maxDpStep_MH $seed_MH $prior $alpha

python3 BayesWHAM_plotter.py

* **BayesReweight.py**

*Reweighting of the biased simulation data and the free energy shifts calculated from solution of the generalized WHAM equations under the Dirichlet prior to estimate F(phi,psi). Phi and psi are backbone dihedral angles with range [-180 deg, 180 deg] and 360 deg periodicity. Within this projection we anticipate good sampling in the driven variable phi but not so in psi where exploration is reliant on spontaneous thermal fluctuations.*

*Copies of all files generated by this calculation are provided in ./diala_phi_EXAMPLE/BayesReweight_OUTPUT*

***Execution:***

T=298;

dim_UMB=1;

periodicity_UMB=[1];

harmonicBiasesFile_UMB='./diala_phi_EXAMPLE/bias/harmonic_biases.txt';

trajDir_UMB='./diala_phi_EXAMPLE/traj';

histBinEdgesFile_UMB='./diala_phi_EXAMPLE/hist/hist_binEdges.txt';

histDir_UMB='./diala_phi_EXAMPLE/hist';

fMAPFile_UMB='f_MAP.txt';

fMHFile_UMB='f_MH.txt';

trajDir_PROJ='./diala_phi_EXAMPLE/traj_aux__phi_psi';

histBinEdgesFile_PROJ='./diala_phi_EXAMPLE/hist_binEdges_phi_psi.txt';

python3 BayesReweight.py $T $dim_UMB $periodicity_UMB $harmonicBiasesFile_UMB $trajDir_UMB $histBinEdgesFile_UMB $histDir_UMB $fMAPFile_UMB $fMHFile_UMB $trajDir_PROJ $histBinEdgesFile_PROJ

python3 BayesReweight_plotter.py



### Example 2 - alanine dipeptide 2D umbrella sampling in [phi,psi] dihedrals###

* **Umbrella sampling input data**

*Ensemble of 36 2D umbrella sampling calculations in [phi,psi] dihedrals*

./diala_phi_psi_EXAMPLE/bias/harmonic_biases.txt - list of umbrella simulation ID, location of harmonic biasing centers (degrees), and harmonic force constants (kJ/mol.degrees^2), one entry per line

./diala_phi_psi_EXAMPLE/traj - trajectories of umbrella sampling variables [phi,psi] in each umbrella sampling calculation

./diala_phi_psi_EXAMPLE/hist/hist_binEdges.txt - bin edges of rectilinear grid used to create [phi,psi] histograms over each umbrella sampling trajectory

./diala_phi_psi_EXAMPLE/hist - [phi,psi] histograms compiled over [phi,psi] trajectories in ./diala_phi_psi_EXAMPLE/traj in row major order

./diala_phi_psi_EXAMPLE/traj_aux__phi_psi_theta - trajectories of projection variables [phi,psi,theta] into which we wish to reweight the FES recorded over the course of each umbrella sampling calculation at the same frequency as those in ./diala_phi_psi_EXAMPLE/traj; the projection variables are arbitrary and need not contain the umbrella sampling variables

* **BayesWHAM.py**

*Calculation of 2D FES F(phi,psi) for alanine dipeptide from 36 2D umbrella sampling simulations at values of phi=[-170:20:170] and psi=[150:20:130]. Phi and psi are backbone dihedral angles with range [-180 deg, 180 deg] and 360 deg periodicity.* 

*We solve the generalized WHAM equations to a tolerance of 1E-15 under a uniform prior. Uncertainties are estimated by performing 1E6 rounds of Metropolis-Hastings sampling from the Bayes posterior and saving every 1E3 realizations. In practice, many more rounds of sampling should be performed until the log likelihood stabilizes indicating the Markov chain burn-in period has passed and samples are being drawn from the stationary distribution.*

*Copies of all files generated by this calculation are provided in ./diala_phi_psi_EXAMPLE/BayesWHAM_OUTPUT*

***Execution:***

dim=2;

periodicity=[1,1];

T=298;

harmonicBiasesFile='./diala_phi_psi_EXAMPLE/bias/harmonic_biases.txt';

histBinEdgesFile='./diala_phi_psi_EXAMPLE/hist/hist_binEdges.txt';

histDir='./diala_phi_psi_EXAMPLE/hist';

tol_WHAM=1E-15;

maxIter_WHAM=1E6;

steps_MH=1E6;

saveMod_MH=1E3;

printMod_MH=1E3;

maxDpStep_MH=5E-5;

seed_MH=200184;

prior='none';

alpha=0;

python3 BayesWHAM.py $dim $periodicity $T $harmonicBiasesFile $histBinEdgesFile $histDir $tol_WHAM $maxIter_WHAM $steps_MH $saveMod_MH $printMod_MH $maxDpStep_MH $seed_MH $prior $alpha

python3 BayesWHAM_plotter.py


* **BayesReweight.py**

*Reweighting of the biased simulation data and the free energy shifts calculated from solution of the generalized WHAM equations under the uniform prior to estimate F(phi,psi,theta). Phi, psi, and theta are backbone dihedral angles with range [-180 deg, 180 deg] and 360 deg periodicity. Within this projection we anticipate good sampling in the driven variables phi and psi but not so in theta where exploration is reliant on spontaneous thermal fluctuations.*

*Copies of all files generated by this calculation are provided in ./diala_phi_psi_EXAMPLE/BayesReweight_OUTPUT*

***Execution:***

T=298;

dim_UMB=2;

periodicity_UMB=[1,1];

harmonicBiasesFile_UMB='./diala_phi_psi_EXAMPLE/bias/harmonic_biases.txt';

trajDir_UMB='./diala_phi_psi_EXAMPLE/traj';

histBinEdgesFile_UMB='./diala_phi_psi_EXAMPLE/hist/hist_binEdges.txt';

histDir_UMB='./diala_phi_psi_EXAMPLE/hist';

fMAPFile_UMB='f_MAP.txt';

fMHFile_UMB='f_MH.txt';

trajDir_PROJ='./diala_phi_psi_EXAMPLE/traj_aux__phi_psi_theta';

histBinEdgesFile_PROJ='./diala_phi_psi_EXAMPLE/hist_binEdges_phi_psi_theta.txt';

python3 BayesReweight.py $T $dim_UMB $periodicity_UMB $harmonicBiasesFile_UMB $trajDir_UMB $histBinEdgesFile_UMB $histDir_UMB $fMAPFile_UMB $fMHFile_UMB $trajDir_PROJ $histBinEdgesFile_PROJ

python3 BayesReweight_plotter.py

### Citing ###

If you use the BayesWHAM.py and BayesReweight.py codes, please read and cite:

A.L. Ferguson "BayesWHAM: A Bayesian approach for free energy estimation, reweighting, and uncertainty quantification in the weighted histogram analysis method" *J. Comput. Chem.* **38** 18 1583-1605 (2017)


### Contact ###

Andrew Ferguson, PhD

Assistant Professor of Materials Science and Engineering

University of Illinois at Urbana-Champaign


[alf@illinois.edu](mailto:alf@illinois.edu)

[http://ferguson.matse.illinois.edu](http://ferguson.matse.illinois.edu)
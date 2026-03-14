"""

% Copyright:	Andrew L. Ferguson, UIUC 
% Last updated:	2 Jan 2016

% SYNOPSIS
%
% code to compute dim-dimensional free energy surface from umbrella simulations by maximization of Bayes posterior, where posterior = likelihood * prior / evidence ~ likelihood * pior, plus estimate uncertainties by Metropolis-Hastings (MH) sampling from posterior  
%
% 1. assumes 	(i)   harmonic restraining potentials in dim-dimensional umbrella variables psi 
%            	(ii)  all simulations conducted at same temperature, T
%            	(iii) rectilinear binning of histograms collected over each umbrella simulation  
% 2. computes maximum a posteriori (MAP) estimate of unbiased probability distribution p(psi) by maximizing P(theta|data) = P(data|theta)*P(theta)/P(data), which is equivalent to solving WHAM equation with a Bayesian prior; 
%    with no prior (i.e., uniform prior) this is the maximum likelihood estimate, which is precisely equivalent to solving the WHAM equations  
% 3. estimates uncertainties by sampling from the posterior distribution P(theta|data) using Metropolis-Hastings algorithm 

% REFERENCES
%
% WHAM:				Kumar, S. et al., J. Comput. Chem. 13 8 1011-1021 (1992) 
%					Kumar, S. et al., J. Comput. Chem. 16 11 1339-1350 (1995) 
%					Ferrenberg, A.M. & Swendsen, R.H., Phys. Rev. Lett. 63 1195-1198 (1989) 
%					Roux, B., Comput. Phys. Comm. 91 275-282 (1995) 
% Bayes sampling:	Gallichio, E. et al., J. Phys. Chem. B 109 6722-6731 (2005)	
%					Bartels, C., Chem. Phys. Lett. 331 446-454 (2000) 
%					Zhu, F. & Hummer, G., J. Comput. Chem. 33 4 453-465 (2012) 
%					Habeck, M., Phys. Rev. Lett. 109 100601 (2012) 

% INPUTS
%
% dim                       - [int] dimensionality of umbrella sampling data in psi(1:dim) = number of coordinates in which umbrella sampling was conducted 
% periodicity               - [1 x dim bool] periodicity in each dimension across range of histogram bins 
%                             -> periodicity(i) == 0 => no periodicity in dimension i; periodicity(i) == 1 => periodicity in dimension i with period defined by range of histogram bins 
% T                         - [float] temperature in Kelvin at which all umbrella sampling simulations were conducted  
% harmonicBiasesFile        - [str] path to text file containing (1 + 2*dim) columns and S rows listing location and strength of biasing potentials in each of the S umbrella simulations 
%                             -> col 1 = umbrella simulation index 1..S, col 2:(2+dim-1) = biased simulation umbrella centers / (arbitrary units), col (2+dim):(2+2*dim-1) = harmonic restraining potential force constant / kJ/mol.(arbitrary units)^2 
% histBinEdgesFile          - [str] path to text file containing dim rows specifying edges of the rectilinear histogram bins in each dimension used to construct dim-dimensional histograms held in histDir/hist_i.txt 
%							  -> histBinEdgesFile contains k=1..dim lines each holding a row vector specifying the rectilinear histogram bin edges in dimension k 
%                             -> for M_k bins in dimension k, there are (M_k+1) edges 
%                             -> bins need not be regularly spaced
% histDir                   - [str] path to directory holding i=1..S dim-dimensional histograms in files hist_i.txt compiled from S biased trajectories over dim-dimensional rectilinear histogram grid specified in histBinEdgesFile  
%                             -> hist_i.txt comprises a row vector containing product_k=1..dim M_k = (M_1*M_2*...*M_dim) values recording histogram counts in each bin of the rectilinear histogram 
%                             -> values recorded in row major order (last index changes fastest) 
% tol_WHAM                  - [float] tolerance on maximum change in any element of p(psi) in self-consistent WHAM solution / - 
% maxIter_WHAM              - [int] maximum # steps allowed to converge to MAP probability estimate by direct iteration of WHAM equations 
% steps_MH                  - [int] # steps in Metropolis-Hastings sampling of Bayes posterior in p(psi) 
% saveMod_MH                - [int] Metropolis-Hastings save modulus 
% printMod_MH               - [int] Metropolis-Hastings print to screen modulus 
% maxDpStep_MH              - [float] maximum size of uniformly distributed step in p(psi) in randomly selected bin 
%                             -> tune to get acceptance probabilities of ~5-20% is a good rule of thumb 
% seed_MH					- [int] seed to initialize random number generator used to perform Metropolis-Hastings sampling of Bayes posterior
% prior                     - [str] type of prior to employ in Bayesian estimate of posterior probability distribution given umbrella sampling histograms
%                             -> prior = {'none','Dirichlet','Gaussian'} 
%                             -> N.B. WHAM equations are solved and prior applied to only the non-zero histogram bins (i.e., in case of Dirichlet pseudo-counts are only added to non-zero count histogram bins) 
% alpha                     - [float] hyperparameter of prior distribution
%                             -> prior = 'Dirichlet' : concentration parmeter for (symmetric) Dirichlet prior, equivalent to adding pseudo-counts of (alpha-1) to each histogram bin; alpha = 1 => uniform prior; alpha > 1 => favors dense + even posteriors, 0 < alpha < 1 => sparse + concentrated posteriors 
%                                        'Gaussian'  : scale parameter for Gaussian penalization of probabilities, exp(-alpha*p^2) 

% OUTPUTS
%
% hist_binCenters.txt       - [dim x M_k float] text file containing dim rows specifying centers of the rectilinear histogram bins in each dimension used to construct dim-dimensional histograms held in histDir/hist_i.txt 
%                             -> dimension k contains M_k bin centers 
%                             -> bins need not be regularly spaced
% hist_binWidths.txt        - [dim x M_k float] text file containing dim rows specifying widths of the rectilinear histogram bins in each dimension used to construct dim-dimensional histograms held in histDir/hist_i.txt 
%                             -> dimension k contains M_k bin widths 
%                             -> bins need not be regularly spaced
% p_MAP.txt                 - [1 x M float] text file containing MAP estimate of unbiased probability distribution p_l over l=1..M bins of dim-dimensional rectilinear histogram grid specified in histBinEdgesFile  
%                             -> values correspond to total probability residing within bin 
%                             -> values recorded in row major order (last index changes fastest) 
% pdf_MAP.txt               - [1 x M float] text file containing MAP estimate of unbiased probability density function pdf_l over l=1..M bins of dim-dimensional rectilinear histogram grid specified in histBinEdgesFile  
%                             -> values correspond to average probability density over the bin 
%                             -> values recorded in row major order (last index changes fastest) 
% betaF_MAP.txt             - [1 x M float] text file containing MAP estimate of unbiased free energy surface beta*F_l = -ln(p(psi)/binVolume) + const. over l=1..M bins of dim-dimensional rectilinear histogram grid specified in histBinEdgesFile to within an arbitrary additive constant 
%                             -> specified only up to an additive constant defined by shifting landscape to have zero mean; when plotting multiple free energy surfaces this is equivalent to minimizing the mean squared error between the landscapes 
%                             -> values recorded in row major order (last index changes fastest) 
% f_MAP.txt                 - [1 x S float] text file containing MAP estimates of f_i = Z/Z_i = ratio of unbiased partition function to that of biased simulation i, for i=1..S biased simulations 
% logL_MAP.txt              - [float] text file containing log of Bayes posterior up to an additive constant of the MAP estimate 
% p_MH.txt                  - [nSamples_MH x M float] text file containing nSamples_MH Metropolis-Hastings samples from the Bayes posterior of unbiased probability distribution p_l over l=1..M bins of dim-dimensional rectilinear histogram grid specified in histBinEdgesFile  
%                             -> values correspond to total probability residing within bin, NOT probability density
%                             -> values recorded in row major order (last index changes fastest) 
% pdf_MH.txt                - [nSamples_MH x M float] text file containing nSamples_MH Metropolis-Hastings samples from the Bayes posterior of unbiased probability density function pdf_l over l=1..M bins of dim-dimensional rectilinear histogram grid specified in histBinEdgesFile  
%                             -> values correspond to average probability density over the bin 
%                             -> values recorded in row major order (last index changes fastest) 
% betaF_MH.txt              - [nSamples_MH x M float] text file containing nSamples_MH Metropolis-Hastings samples from the Bayes posterior of unbiased free energy surface beta*F_l = -ln(p(psi)/binVolume) + const. over l=1..M bins of dim-dimensional rectilinear histogram grid specified in histBinEdgesFile to within an arbitrary additive constant 
%                             -> specified only up to an additive constant defined by shifting landscape to have zero mean; when plotting multiple free energy surfaces this is equivalent to minimizing the mean squared error between the landscapes 
%                             -> values recorded in row major order (last index changes fastest) 
% f_MH.txt                  - [nSamples_MH x S float] text file containing nSamples_MH Metropolis-Hastings samples from the Bayes posterior of f_i = Z/Z_i = ratio of unbiased partition function to that of biased simulation i, for i=1..S biased simulations 
% logL_MH.txt               - [nSamples_MH x 1 float] text file containing log of Bayes posterior up to an additive constant over the nSamples_MH Metropolis-Hastings samples from the posterior 
% step_MH.txt               - [nSamples_MH x 1 int] text file containing Metropolis-Hastings step associated with each of the nSamples_MH Metropolis-Hastings samples from the posterior 

"""

## imports
import os, re, sys, time
import random, math
import yaml

import numpy as np
import numpy.matlib
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

# Python version check
if sys.version_info[0] < 3:
    raise Exception("Python 3 or a more recent version is required.")

## classes

## methods

# usage
def _usage():
	print("USAGE: %s [--config CONFIG_FILE | LEGACY_ARGS...]" % sys.argv[0])
	print("")
	print("  YAML Configuration (Recommended):")
	print("       --config CONFIG_FILE      - Path to YAML configuration file")
	print("                                   (see config_example.yaml for template)")
	print("")
	print("  Legacy Command-Line Arguments (for backward compatibility):")
	print("       %s dim periodicity T harmonicBiasesFile histBinEdgesFile histDir tol_WHAM maxIter_WHAM steps_MH saveMod_MH printMod_MH maxDpStep_MH seed_MH prior alpha [kB]" % sys.argv[0])
	print("")
	print("       dim                       - [int] dimensionality of umbrella sampling data in psi(1:dim) = number of coordinates in which umbrella sampling was conducted ")
	print("       periodicity               - [1 x dim bool] periodicity in each dimension across range of histogram bins ")
	print("                                   -> periodicity(i) == 0 => no periodicity in dimension i; periodicity(i) == 1 => periodicity in dimension i with period defined by range of histogram bins ")
	print("       T                         - [float] temperature in Kelvin at which all umbrella sampling simulations were conducted ")
	print("       harmonicBiasesFile        - [str] path to text file containing (1 + 2*dim) columns and S rows listing location and strength of biasing potentials in each of the S umbrella simulations ")
	print("                                   -> col 1 = umbrella simulation index 1..S, col 2:(2+dim-1) = biased simulation umbrella centers / (arbitrary units), col (2+dim):(2+2*dim-1) = harmonic restraining potential force constant / kJ/mol.(arbitrary units)^2 ")
	print("       histBinEdgesFile          - [str] path to text file containing dim rows specifying edges of the rectilinear histogram bins in each dimension used to construct dim-dimensional histograms held in histDir/hist_i.txt ")
	print("                                   -> histBinEdgesFile contains k=1..dim lines each holding a row vector specifying the rectilinear histogram bin edges in dimension k ")
	print("                                   -> for M_k bins in dimension k, there are (M_k+1) edges ")
	print("                                   -> bins need not be regularly spaced ")
	print("       histDir                   - [str] path to directory holding i=1..S dim-dimensional histograms in files hist_i.txt compiled from S biased trajectories over dim-dimensional rectilinear histogram grid specified in histBinEdgesFile ")
	print("                                   -> hist_i.txt comprises a row vector containing product_k=1..dim M_k = (M_1*M_2*...*M_dim) values recording histogram counts in each bin of the rectilinear histogram ")
	print("                                   -> values recorded in row major order (last index changes fastest) ")
	print("       tol_WHAM                  - [float] tolerance on maximum change in any element of p(psi) in self-consistent WHAM solution / - ")
	print("       maxIter_WHAM              - [int] maximum # steps allowed to converge to MAP probability estimate by direct iteration of WHAM equations ")
	print("       steps_MH                  - [int] # steps in Metropolis-Hastings sampling of Bayes posterior in p(psi) ")
	print("       saveMod_MH                - [int] Metropolis-Hastings save modulus ")
	print("       printMod_MH               - [int] Metropolis-Hastings print to screen modulus ")
	print("       maxDpStep_MH              - [float] maximum size of uniformly distributed step in p(psi) in randomly selected bin ")
	print("                                   -> tune to get acceptance probabilities of ~5-20% is a good rule of thumb ")
	print("       seed_MH					- [int] seed to initialize random number generator used to perform Metropolis-Hastings sampling of Bayes posterior ")
	print("       prior                     - [str] type of prior to employ in Bayesian estimate of posterior probability distribution given umbrella sampling histograms ")
	print("                                   -> prior = {'none','Dirichlet','Gaussian'} ")
	print("                                   -> N.B. WHAM equations are solved and prior applied to only the non-zero histogram bins (i.e., in case of Dirichlet pseudo-counts are only added to non-zero count histogram bins) ")
	print("       alpha                     - [float] hyperparameter of prior distribution ")
	print("                                   -> prior = 'Dirichlet' : concentration parmeter for (symmetric) Dirichlet prior, equivalent to adding pseudo-counts of (alpha-1) to each histogram bin; alpha = 1 => uniform prior; alpha > 1 => favors dense + even posteriors, 0 < alpha < 1 => sparse + concentrated posteriors ")
	print("                                              'Gaussian'  : scale parameter for Gaussian penalization of probabilities, exp(-alpha*p^2) ")
	print("       kB (optional)             - [float] Boltzmann constant in desired energy units (default: 0.0083144621 kJ/mol·K)")
	print("                                   -> Common values: 0.0083144621 (kJ/mol·K), 0.001987204 (kcal/mol·K), 8.314462618e-3 (kJ/mol·K), 1.380649e-23 (J/K)")

def load_config_from_yaml(config_file):
	"""Load configuration from YAML file"""
	with open(config_file, 'r') as f:
		config = yaml.safe_load(f)

	# Extract parameters from YAML structure
	params = {}
	params['dim'] = config['umbrella_sampling']['dim']
	params['periodicity'] = config['umbrella_sampling']['periodicity']
	params['T'] = config['umbrella_sampling']['temperature']
	params['kB'] = config['umbrella_sampling'].get('boltzmann_constant', 0.0083144621)  # Default: kJ/mol·K
	params['harmonicBiasesFile'] = config['input_files']['harmonic_biases']
	params['histBinEdgesFile'] = config['input_files']['hist_bin_edges']
	params['histDir'] = config['input_files']['hist_dir']
	params['tol_WHAM'] = config['wham']['tolerance']
	params['maxIter_WHAM'] = config['wham']['max_iterations']
	params['steps_MH'] = config['metropolis_hastings']['steps']
	params['saveMod_MH'] = config['metropolis_hastings']['save_modulus']
	params['printMod_MH'] = config['metropolis_hastings']['print_modulus']
	params['maxDpStep_MH'] = config['metropolis_hastings']['max_step_size']
	params['seed_MH'] = config['metropolis_hastings']['random_seed']
	params['prior'] = config['prior']['type']
	params['alpha'] = config['prior']['alpha']
	params['output_dir'] = config.get('output_dir', None)

	return params	
	
def unique(a):
    """ return the list with duplicate elements removed """
    return list(set(a))

def intersect(a, b):
    """ return the intersection of two lists """
    return list(set(a) & set(b))

def union(a, b):
    """ return the union of two lists """
    return list(set(a) | set(b))

def ind2sub_RMO(sz,ind):

	if (ind < 0) or (ind > (np.prod(sz)-1)):
		print("\nERROR - Variable ind = %d < 0 or ind = %d > # elements in tensor of size sz = %d in ind2sub_RMO" % (ind,ind,np.prod(sz)))
		sys.exit(-1)
	
	sub = np.zeros(sz.shape[0], dtype=np.uint64)
	for ii in range(0,len(sz)-1):
		P = np.prod(sz[ii+1:])
		sub_ii = math.floor(float(ind)/float(P))
		sub[ii] = sub_ii
		ind = ind - sub_ii*P
	sub[-1] = ind
	
	return sub

def sub2ind_RMO(sz,sub):
	
	if sz.shape != sub.shape:
		print("\nERROR - Variables sub and sz in sub2ind_RMO are not commensurate")
		sys.exit(-1)
	
	for ii in range(0,len(sz)-1):
		if ( sub[ii] < 0 ) or ( sub[ii] >= sz[ii] ):
			print("\nERROR - sub[%d] = %d < 0 or sub[%d] = %d >= sz[%d] = %d in sub2ind_RMO"  % (ii,sub[ii],ii,sub[ii],ii,sz[ii]))
			sys.exit(-1)
	
	ind=0
	for ii in range(0,len(sz)-1):
		ind = ind + sub[ii]*np.prod(sz[ii+1:])
	ind = ind + sub[-1]

	return int(ind)


## main

# parameters
kB = 0.0083144621		# Boltzmann's constant / kJ/mol.K
genGaussPriorPower = 2	# power of generalized Gaussian implemented for prior='Gaussian' (genGaussPriorPower=1 => Laplacian, genGaussPriorPower=2 => Gaussian, etc.)
						# - hard coded as standard Gaussian prior with genGaussPriorPower = 2, but other values are possible
						# - Laplacian prior has no impact (i.e., same as uniform prior) due to normalization of probability distribution since appears in log posterior as nothing more than an additive constant = -alpha  

# loading inputs
# Check if using YAML config or legacy CLI arguments
if len(sys.argv) == 3 and sys.argv[1] == '--config':
	# YAML configuration mode
	config_file = sys.argv[2]
	print("Loading configuration from: %s" % config_file)
	params = load_config_from_yaml(config_file)

	dim = params['dim']
	periodicity = params['periodicity']
	T = params['T']
	kB = params['kB']
	harmonicBiasesFile = params['harmonicBiasesFile']
	histBinEdgesFile = params['histBinEdgesFile']
	histDir = params['histDir']
	tol_WHAM = params['tol_WHAM']
	maxIter_WHAM = params['maxIter_WHAM']
	steps_MH = params['steps_MH']
	saveMod_MH = params['saveMod_MH']
	printMod_MH = params['printMod_MH']
	maxDpStep_MH = params['maxDpStep_MH']
	seed_MH = params['seed_MH']
	prior = params['prior']
	alpha = params['alpha']
	output_dir = params['output_dir']

	# Change to output directory if specified
	if output_dir is not None:
		os.makedirs(output_dir, exist_ok=True)
		os.chdir(output_dir)
		print("Output will be written to: %s\n" % output_dir)

elif len(sys.argv) == 16 or len(sys.argv) == 17:
	# Legacy CLI arguments mode (15 or 16 args)
	dim = int(sys.argv[1])
	periodicity = str(sys.argv[2])
	T = float(sys.argv[3])
	harmonicBiasesFile = str(sys.argv[4])
	histBinEdgesFile = str(sys.argv[5])
	histDir = str(sys.argv[6])
	tol_WHAM = float(sys.argv[7])
	maxIter_WHAM = int(float(sys.argv[8]))
	steps_MH = int(float(sys.argv[9]))
	saveMod_MH = int(float(sys.argv[10]))
	printMod_MH = int(float(sys.argv[11]))
	maxDpStep_MH = float(sys.argv[12])
	seed_MH = int(sys.argv[13])
	prior = str(sys.argv[14])
	alpha = float(sys.argv[15])
	# Optional kB parameter (16th argument)
	if len(sys.argv) == 17:
		kB = float(sys.argv[16])
	else:
		kB = 0.0083144621  # Default: kJ/mol·K
else:
	_usage()
	sys.exit(-1)

# - post-processing and error checking
beta = 1/(kB*T)			# beta = 1/kBT = (kJ/mol)^(-1)

# Handle periodicity (already a list if from YAML, string if from CLI)
if isinstance(periodicity, str):
	periodicity = periodicity[1:-1].split(',')
	periodicity = [int(x) for x in periodicity]

if prior not in ['none','Dirichlet','Gaussian']:
	print("\nERROR - Input prior = '%s' unrecognized, must be one of {'none','Dirichlet','Gaussian'}" % (prior))
	sys.exit(-1)

# - printing args to screen
print("")
print("dim = %d" % (dim))
print("periodicity = %s" % (periodicity))
print("T = %e" % (T))
print("kB = %e" % (kB))
print("harmonicBiasesFile = %s" % (harmonicBiasesFile))
print("histBinEdgesFile = %s" % (histBinEdgesFile))
print("histDir = %s" % (histDir))
print("tol_WHAM = %e" % (tol_WHAM))
print("maxIter_WHAM = %d" % (maxIter_WHAM))
print("steps_MH = %d" % (steps_MH))
print("saveMod_MH = %d" % (saveMod_MH))
print("printMod_MH = %d" % (printMod_MH))
print("maxDpStep_MH = %e" % (maxDpStep_MH))
print("seed_MH = %d" % (seed_MH))
print("prior = %s" % (prior))
print("alpha = %e" % (alpha))
print("")


# loading data
print("Loading data...")

# - loading location and strength of harmonic biasing potentials in each of S umbrella simulations 
harmonicBiasesFile_DATA = []
with open(harmonicBiasesFile,'r') as fin:
	for line in fin:
		harmonicBiasesFile_DATA.append(line.strip().split())

S = len(harmonicBiasesFile_DATA)

if len(harmonicBiasesFile_DATA[0]) != (1+2*dim):
	print("\nERROR - Number of columns in %s != 1+2*dim = %d as expected" % (harmonicBiasesFile,1+2*dim))
	sys.exit(-1)
	
umbC = [item[1:1+dim] for item in harmonicBiasesFile_DATA]
umbC = [[float(y) for y in x] for x in umbC]
umbC = np.array(umbC)

umbF = [item[1+dim:1+2*dim] for item in harmonicBiasesFile_DATA]
umbF = [[float(y) for y in x] for x in umbF]
umbF = np.array(umbF)

# - loading histogram bin edges 
binE = []
with open(histBinEdgesFile,'r') as fin:
	for line in fin:
		binE.append(line.strip().split())

if len(binE) != dim:
	print("\nERROR - %s does not contain expected number of dim = %d lines specifying histogram bin edges" % (histBinEdgesFile,dim))
	sys.exit(-1)

binE = [[float(y) for y in x] for x in binE]
binE = np.array(binE)

# - counting number of bins in each dimension M_k k=1..dim, and total number of bins M = product_k=1..dim M_k = M_1*M_2*...*M_dim 
M_k = np.zeros(dim, dtype=np.uint64)
for d in range(0,dim):
    M_k[d] = len(binE[d])-1
M = np.prod(M_k)

# - converting binEdges into binCenters and binWidths
binC = []
binW = []
for d in range(0,dim):
	binC_d = np.zeros(M_k[d])
	binW_d = np.zeros(M_k[d])
	for ii in range(0,M_k[d]):
		binC_d[ii] = 0.5 * ( binE[d][ii] + binE[d][ii+1] )
		binW_d[ii] = binE[d][ii+1] - binE[d][ii]
	binC.append(binC_d)
	binW.append(binW_d)

# - computing period in each periodic dimension as histogram bin range 
period = np.zeros(dim)
for d in range(0,dim):
    if periodicity[d] == 0:
    	period[d] = float('nan')
    else:
        period[d] = binE[d][-1] - binE[d][0]


# loading histograms compiled over the S umbrella simulations recorded as row vectors in row major order (last index changes fastest) 
n_il = []
for i in range(0,S):
	hist_filename = histDir + '/hist_' + str(i+1) + '.txt'
	hist_DATA= []
	with open(hist_filename,'r') as fin:
		for line in fin:
			hist_DATA.append(line.strip().split())
	if len(hist_DATA) != 1:
		print("\nERROR - Did not find expected row vector in reading %s" % (hist_filename))
		sys.exit(-1)
	if len(hist_DATA[0]) != M:
		print("\nERROR - Row vector in %s did not contain expected number of elements M = M_1*M_2*...*M_dim = %d given histogram bins specified in %s" % (hist_filename,M,histBinEdgesFile))
		sys.exit(-1)
	n_il.append(hist_DATA[0])
n_il = [[float(y) for y in x] for x in n_il]
n_il = np.array(n_il)


# precomputing aggregated statistics
N_i = np.sum(n_il,axis=1)         # total counts in simulation i
M_l = np.sum(n_il,axis=0)         # total counts in bin l

print("DONE!\n\n")


# checking M_l for gaps 
# -> gaps in the aggregated M_l histogram (i.e., empty bins between populated bins) can cause numerical stability problems for the WHAM equations 
#    - conceptually, the presence of gaps indicates that there is not mutual (at least pairwise) overlap between all biased histograms, meaning that the estimate of the unbiased distribution may be disjoint since the biased distributions cannot be stitched together into a continuous probability distribution 
#    - mathematically, p_l becomes disjoint and cannot be consistently normalized across the full domain resulting in multiple disjoint distributions that leads to numerical instability in the WHAM equations 
# -> contiguity of the histogram (i.e., the absence of empty bins between populated bins) is a sufficient, but not necessary, condition to assure stability of the WHAM equations since this guarantees mutual overlap between all of the biased histograms and a non-disjoint estimate of the unbiased probability distribution (sufficiently small gaps are tolerable provided that the probability distribution does not become disjoint within machine precision) 

# - computing adjacency matrix representing the graph specifying contiguity of populated histogram bins 
'''
# TEST BLOCK for adjacency routines
dim = 3
periodicity = [0,0,0]

M_k = np.array([3,4,2], dtype=np.uint64)
M = np.prod(M_k)

A = [ [ [1,1], [0,0], [1,1], [1,1] ], [ [0,0], [0,0], [0,0], [1,1] ], [ [1,1], [0,0], [1,1], [1,1] ] ]
A=np.array(A)
print A

M_l=np.zeros(M)
for i in range(0,M_k[0]):
	for j in range(0,M_k[1]):
		for k in range(0,M_k[2]):
			#print("%d %d %d\n" % (i,j,k))
			M_l[i*M_k[1]*M_k[2] + j*M_k[2] + k] = A[i,j,k]
'''

print("Checking connectivity of aggregated biased histogram...")
adjacency = np.zeros((M,M))
for l in range(0,M):
	if M_l[l] > 0:
		sub = ind2sub_RMO(M_k,l)
		for d in range(0,dim):
			for q in [-1,+1]:
				
				sub_d_plus_q = sub[d] + q
				if sub_d_plus_q == M_k[d]:
					if periodicity[d] == 0:
						continue
					else:
						sub_d_plus_q = 0			# wrapping through periodic boundary
				elif sub_d_plus_q == -1:
					if periodicity[d] == 0:
						continue
					else:
						sub_d_plus_q = M_k[d]-1		# wrapping through periodic boundary
				sub_d_plus_q = int(sub_d_plus_q)
				
				sub_neigh = np.array(sub)
				sub_neigh[d] = sub_d_plus_q
				
				ind_neigh = sub2ind_RMO(M_k,sub_neigh)
				if M_l[ind_neigh] > 0:
					adjacency[l,ind_neigh] = 1
				
# - applying Tarjan's algorithm to identify contiguous regions of non-zero histogram bins within the adjacency matrix 
maskNZ = (M_l>0)

adjacency_NZ = adjacency[:,maskNZ]
adjacency_NZ = np.transpose(np.transpose(adjacency_NZ)[:,maskNZ])

adjacency_NZ_SPARSE = csr_matrix(adjacency_NZ)

[nConComp,conComp_labels] = connected_components(adjacency_NZ_SPARSE, directed=False, connection='weak', return_labels=True)

if nConComp > 1:
    print("\nWARNING - aggregated biased histogram is disconnected, containing %d non-contiguous connected components\n" % (nConComp))
    print("        - we might anticipate WHAM convergence difficulties\n")
    print("        - consider making the aggregated histogram contiguous by increasing the bin size or conducting additional biased runs to patch the gaps\n")
    time.sleep(1)
print("DONE!\n\n")


# precomputing biasing potentials for each simulation i in each bin l 
print("Computing biasing potentials for each simulation in each bin...")
c_il = np.zeros((S,M))		# S-by-M matrix containing biases due to artificial harmonic potential for simulation i in bin l
for i in range(0,S):
	for l in range(0,M):
		sub = ind2sub_RMO(M_k,l)
		expArg = 0
		for d in range(0,dim):
			if periodicity[d] != 0:
				delta = min( abs(binC[d][sub[d]]-umbC[i,d]), abs(binC[d][sub[d]]-umbC[i,d]+period[d]), abs(binC[d][sub[d]]-umbC[i,d]-period[d]) )
			else:
				delta = abs(binC[d][sub[d]]-umbC[i,d])
			expArg = expArg + 0.5*umbF[i,d]*math.pow(delta,2)
		c_il[i,l] = math.exp(-beta*expArg)
	print("\tProcessing of simulation %d of %d complete" % (i+1,S))
print("DONE!\n\n")

print("\n=== BIAS POTENTIAL MATRIX DIAGNOSTICS ===")
print("c_il matrix shape: (%d, %d)" % c_il.shape)
print("c_il statistics: min=%e, max=%e, mean=%e" % (c_il.min(), c_il.max(), c_il.mean()))
print("Checking for problematic simulations:")
problematic_sims = []
for i in range(S):
	row_nonzero = np.sum(c_il[i,:] > 0)
	if row_nonzero < 10:
		print("  WARNING: Sim %d has only %d/%d non-zero bias values" % (i+1, row_nonzero, M))
		problematic_sims.append(i+1)
	elif i < 3:  # Always show first 3
		print("  Sim %d: %d/%d non-zero bias values" % (i+1, row_nonzero, M))
if len(problematic_sims) > 0:
	print("  Total problematic simulations: %s" % str(problematic_sims))
zero_bias_bins = np.sum(np.all(c_il == 0, axis=0))
if zero_bias_bins > 0:
	print("WARNING: %d bins have zero bias from all simulations" % zero_bias_bins)

# MAP estimate of unbiased probability distribution over the umbrella histogram bins p_l by direct iteration of the WHAM equations
# -> analogy of WHAM equations under Bayesian prior derived using method of Lagrangian multipliers 
# -> Ref: C. Bartels "Analyzing biased Monte Carlo and molecular dynamics simulations" Chemical Physics Letters 331 446-454 (2000)

# - masking out all zero-count bins not visited by any simulation for computational efficiency 
#   -> optimal probability estimate for unvisited bins with M_l=0 is p_l=0, so can eliminate these bins from consideration 
#   -> Ref: C. Bartels "Analyzing biased Monte Carlo and molecular dynamics simulations" Chemical Physics Letters 331 446-454 (2000)
maskNZ = (M_l>0)
M_l_NZ = np.array(M_l[maskNZ])
c_il_NZ = np.array(c_il[:,maskNZ])
M_NZ = M_l_NZ.size

print("\n=== INITIAL STATE DIAGNOSTICS ===")
print("Number of simulations (S): %d" % S)
print("Total bins (M): %d" % M)
print("Non-zero bins (M_NZ): %d" % M_NZ)
print("N_i (counts per simulation):")
for i in range(S):
	print("  Simulation %d: %d counts" % (i+1, N_i[i]))
print("M_l_NZ statistics:")
print("  Min: %d, Max: %d, Mean: %.2f" % (M_l_NZ.min(), M_l_NZ.max(), M_l_NZ.mean()))
print("c_il_NZ matrix statistics:")
print("  Min: %e, Max: %e" % (c_il_NZ.min(), c_il_NZ.max()))
print("  Rows with all zeros: %d" % np.sum(np.all(c_il_NZ == 0, axis=1)))
print("  Columns with all zeros: %d" % np.sum(np.all(c_il_NZ == 0, axis=0)))

# - initial guess for p_l and f_i
p_l_NZ = np.ones(M_NZ)
p_l_NZ = np.divide( p_l_NZ, np.sum(p_l_NZ) )
f_i = np.divide( 1, np.dot(c_il_NZ,p_l_NZ) )

print("\n=== INITIAL GUESS ===")
print("p_l_NZ: sum=%.10f, min=%e, max=%e" % (np.sum(p_l_NZ), p_l_NZ.min(), p_l_NZ.max()))
dot_product = np.dot(c_il_NZ, p_l_NZ)
print("Dot product c_il_NZ @ p_l_NZ:")
for i in range(S):
	print("  Sim %d: %e" % (i+1, dot_product[i]))
print("Initial f_i values:")
for i in range(S):
	if dot_product[i] == 0:
		print("  Sim %d: INFINITE (zero overlap detected!)" % (i+1))
	else:
		print("  Sim %d: %e" % (i+1, f_i[i]))
zero_overlap_count = np.sum(dot_product == 0)
if zero_overlap_count > 0:
	print("WARNING: %d simulations have zero overlap!" % zero_overlap_count)
print("Simulations with f_i > 1e10: %d" % np.sum(f_i > 1e10))

# - direct iteration of WHAM equations to self-consistency at tolerance = tol 
dp_max = float('Inf')
dlogf_max = float('Inf')
iter = 1
print("Computing MAP estimate of probability distribution by direct iteration of the WHAM equations...")
while (dp_max > tol_WHAM):

	p_l_NZ_OLD = np.array(p_l_NZ)
	f_i_OLD = np.array(f_i)

	if iter <= 5 or iter % 10 == 0:  # Print first 5 iterations and every 10th
		print("\n--- Iteration %d Diagnostics ---" % iter)
		print("p_l_NZ_OLD: min=%e, max=%e, any NaN=%s, any Inf=%s" %
		      (p_l_NZ_OLD.min(), p_l_NZ_OLD.max(),
		       str(np.any(np.isnan(p_l_NZ_OLD))), str(np.any(np.isinf(p_l_NZ_OLD)))))
		print("f_i_OLD: min=%e, max=%e, any NaN=%s, any Inf=%s" %
		      (f_i_OLD.min(), f_i_OLD.max(),
		       str(np.any(np.isnan(f_i_OLD))), str(np.any(np.isinf(f_i_OLD)))))

		denominator_preview = np.dot(np.multiply(N_i, f_i), c_il_NZ)
		print("Denominator for p_l_NZ update:")
		print("  min=%e, max=%e" % (denominator_preview.min(), denominator_preview.max()))
		num_zeros = np.sum(denominator_preview == 0)
		if num_zeros > 0:
			print("  WARNING: %d bins have zero denominator!" % num_zeros)
			zero_indices = np.where(denominator_preview == 0)[0]
			print("  Zero denominator at bin indices: %s" % str(zero_indices[:10]))

	if prior == 'none':
		p_l_NZ = np.divide( M_l_NZ, np.dot( np.multiply(N_i,f_i), c_il_NZ ) )
	elif prior == 'Dirichlet':
		p_l_NZ = np.divide( ( M_l_NZ + (alpha-1)*np.ones(M_NZ) ), np.dot( np.multiply(N_i,f_i), c_il_NZ ) + M_NZ*(alpha-1)*np.ones(M_NZ) )
	elif prior == 'Gaussian':
		if alpha == 0:		# no iterative solution required for alpha=0 => uniform prior => no prior
			p_l_NZ = np.divide( M_l_NZ, np.dot( np.multiply(N_i,f_i), c_il_NZ ) )
		else:
			a = genGaussPriorPower*alpha*np.ones(M_NZ)
			b = np.dot( np.multiply(N_i,f_i), c_il_NZ )
			c = M_l_NZ

			dp_max_inner = float('Inf')
			p_l_NZ_inner = np.array(p_l_NZ)
			while (dp_max_inner > tol_WHAM):		# iterative solution for regularized probability distribution under generalized Gaussian prior to same tolerance as WHAM equations 
				p_l_NZ_inner_OLD = np.array(p_l_NZ_inner)
				p_l_NZ_inner_POWt = np.array([math.pow(x,genGaussPriorPower) for x in p_l_NZ_inner])
				p_l_NZ_inner = np.divide( c - np.multiply( a, p_l_NZ_inner_POWt ), b - np.dot( np.ones(M_NZ), p_l_NZ_inner_POWt )*a )
				dp_max_inner = max( abs(p_l_NZ_inner - p_l_NZ_inner_OLD) )
			p_l_NZ = np.array(p_l_NZ_inner)
	else:
		print("\nERROR - Input prior = '%s' unrecognized, must be one of {'none','Dirichlet','Gaussian'}" % (prior))
		sys.exit(-1)
	p_l_NZ = np.divide( p_l_NZ, np.sum(p_l_NZ) )

	if iter <= 5 or iter % 10 == 0:
		print("Updated p_l_NZ:")
		print("  sum=%.15f, min=%e, max=%e" % (np.sum(p_l_NZ), p_l_NZ.min(), p_l_NZ.max()))
		print("  any NaN=%s, any Inf=%s, num zeros=%d" %
		      (str(np.any(np.isnan(p_l_NZ))), str(np.any(np.isinf(p_l_NZ))), np.sum(p_l_NZ == 0)))

		# Check denominator for f_i update
		f_i_denominator = np.dot(c_il_NZ, p_l_NZ)
		num_zero_denom = np.sum(f_i_denominator == 0)
		if num_zero_denom > 0:
			print("  WARNING: %d simulations will have Inf f_i!" % num_zero_denom)
			zero_sims = np.where(f_i_denominator == 0)[0] + 1
			print("  Problematic simulations: %s" % str(zero_sims))

	f_i = np.divide( 1, np.dot(c_il_NZ,p_l_NZ) )

	if iter <= 5 or iter % 10 == 0:
		print("Updated f_i:")
		print("  min=%e, max=%e" % (f_i.min(), f_i.max()))
		print("  any NaN=%s, any Inf=%s" % (str(np.any(np.isnan(f_i))), str(np.any(np.isinf(f_i)))))
		num_large = np.sum(f_i > 1e10)
		if num_large > 0:
			print("  WARNING: %d simulations have f_i > 1e10" % num_large)
			large_sims = np.where(f_i > 1e10)[0] + 1
			print("  Large f_i simulations: %s" % str(large_sims[:10]))

	dp_max = max( abs(p_l_NZ - p_l_NZ_OLD) )

	print("\n--- Convergence Check (Iteration %d) ---" % iter)
	print("dp_max = %s (is finite: %s)" % (str(dp_max), str(np.isfinite(dp_max))))
	if not np.isfinite(dp_max):
		print("p_l_NZ has %d NaN, %d Inf values" % (np.sum(np.isnan(p_l_NZ)), np.sum(np.isinf(p_l_NZ))))
		print("p_l_NZ_OLD has %d NaN, %d Inf values" % (np.sum(np.isnan(p_l_NZ_OLD)), np.sum(np.isinf(p_l_NZ_OLD))))

	if not np.isfinite(dp_max):
		print("\n" + "="*70)
		print("ERROR - WHAM equation convergence failed at iteration %d" % iter)
		print("="*70)
		print("\nCONVERGENCE METRIC:")
		print("  dp_max = %s" % str(dp_max))
		print("\nPROBABILITY DISTRIBUTION STATE:")
		print("  p_l_NZ: min=%s, max=%s" % (str(p_l_NZ.min()), str(p_l_NZ.max())))
		print("  p_l_NZ has %d NaN values" % np.sum(np.isnan(p_l_NZ)))
		print("  p_l_NZ has %d Inf values" % np.sum(np.isinf(p_l_NZ)))
		print("  p_l_NZ has %d zero values" % np.sum(p_l_NZ == 0))
		print("\nPARTITION FUNCTION RATIOS:")
		print("  f_i: min=%s, max=%s" % (str(f_i.min()), str(f_i.max())))
		print("  f_i has %d NaN values" % np.sum(np.isnan(f_i)))
		print("  f_i has %d Inf values" % np.sum(np.isinf(f_i)))
		if np.any(np.isinf(f_i)):
			inf_sims = np.where(np.isinf(f_i))[0] + 1
			print("  Simulations with infinite f_i: %s" % str(inf_sims))
		print("\nOVERLAP ANALYSIS:")
		overlap_check = np.dot(c_il_NZ, p_l_NZ_OLD)
		zero_overlap = np.sum(overlap_check == 0)
		print("  Simulations with zero overlap: %d of %d" % (zero_overlap, S))
		if zero_overlap > 0:
			zero_sims = np.where(overlap_check == 0)[0] + 1
			print("  Zero overlap simulations: %s" % str(zero_sims))
		print("\nLIKELY CAUSE:")
		if zero_overlap > 0:
			print("  => Insufficient overlap between umbrella windows")
			print("  => Simulations %s have no overlap with sampled regions" % str(zero_sims))
		elif np.any(np.isinf(f_i)):
			print("  => Numerical instability in partition function ratios")
		else:
			print("  => Unknown numerical issue - check input data quality")
		print("\nRECOMMENDATIONS:")
		print("  1. Check umbrella window spacing and overlap")
		print("  2. Verify histogram data quality (no all-zero histograms)")
		print("  3. Consider adjusting force constants or bin sizes")
		print("  4. Review biasing potential configuration")
		print("  5. Try using Dirichlet prior (alpha=2.0) for numerical stability")
		print("="*70 + "\n")
		sys.exit(-1)
	
	log_f_i = np.array([math.log(x) for x in f_i])
	log_f_i_OLD = np.array([math.log(x) for x in f_i_OLD])
	dlogf_max = max( abs( log_f_i - log_f_i_OLD ) )
	
	print("\tIteration %d, max(dp) = %15.5e" % (iter,dp_max))
	
	iter = iter + 1
	if iter > maxIter_WHAM:
		print("\nERROR - WHAM iteration did not converge within maxIter_WHAM = %d" % (maxIter_WHAM))
		sys.exit(-1)
		
print("\tWHAM eqautions converged to within tolerance of %15.5e" % (tol_WHAM))
print("DONE!\n\n")

# - unmasking p_l_NZ into p_l_MAP by reconstituting empty bins l with M_l=0 with optimal values p_l=0 
p_l_MAP = np.zeros(M)
p_l_MAP[maskNZ] = p_l_NZ
f_i_MAP = np.array(f_i)

# - converting unbiased probability distribution estimate p_l into probability density function pdf_l and free energy estimate betaF_l 
#   -> mean zeroing betaF; when plotting multiple free energy surfaces this is equivalent to minimizing the mean squared error between the landscapes 
pdf_l_MAP = np.zeros(p_l_MAP.size)
betaF_l_MAP = np.ones(p_l_MAP.size)*float('nan')
for l in range(0,M):
	sub = ind2sub_RMO(M_k,l)
	binVol = 1
	for d in range(0,dim):
		binVol = binVol*binW[d][sub[d]]
	if p_l_MAP[l] > 0:
		pdf_l_MAP[l] = p_l_MAP[l]/binVol
		betaF_l_MAP[l] = -math.log(p_l_MAP[l]/binVol)
betaF_l_MAP = betaF_l_MAP - np.mean(betaF_l_MAP[np.isfinite(betaF_l_MAP)])

# - writing to file bin centers and bin widths of l=1..M bins of dim-dimensional rectilinear histogram grid as row vectors in row major order (last index changes fastest) 
with open("hist_binCenters.txt","w") as fout:
	for d in range(0,dim):
		for ii in range(0,M_k[d]):
			fout.write("%15.5e" % (binC[d][ii]))
		fout.write("\n")

with open("hist_binWidths.txt","w") as fout:
	for d in range(0,dim):
		for ii in range(0,M_k[d]):
			fout.write("%15.5e" % (binW[d][ii]))
		fout.write("\n")

# - writing to file MAP estimates for p_l, pdf_l, and betaF_l over l=1..M bins of dim-dimensional rectilinear histogram grid as row vectors in row major order (last index changes fastest) 
with open("p_MAP.txt","w") as fout:
	for l in range(0,M):
		fout.write("%15.5e" % (p_l_MAP[l]))
	fout.write("\n")

with open("pdf_MAP.txt","w") as fout:
	for l in range(0,M):
		fout.write("%15.5e" % (pdf_l_MAP[l]))
	fout.write("\n")

with open("betaF_MAP.txt","w") as fout:
	for l in range(0,M):
		fout.write("%15.5e" % (betaF_l_MAP[l]))
	fout.write("\n")

# - writing to file MAP estimates for f_i = Z/Z_i = ratio of unbiased partition function to that of biased simulation i,for i=1..S biased simulations 
with open("f_MAP.txt","w") as fout:
	for i in range(0,S):
		fout.write("%15.5e" % (f_i_MAP[i]))
	fout.write("\n")
	
# - writing to file log of Bayes posterior up to an additive constant for the MAP estimate 
with open("logL_MAP.txt","w") as fout:
	maskNZ = (M_l>0)
	M_l_NZ = np.array(M_l[maskNZ])
	p_l_MAP_NZ = np.array(p_l_MAP[maskNZ])
	M_NZ = M_l_NZ.size
	if prior == 'none':
		logL_MAP = np.dot( N_i, np.array([math.log(x) for x in f_i_MAP]) ) + np.dot( M_l_NZ, np.array([math.log(x) for x in p_l_MAP_NZ]) )
	elif prior == 'Dirichlet':
		logL_MAP = np.dot( N_i, np.array([math.log(x) for x in f_i_MAP]) ) + np.dot( M_l_NZ + (alpha-1)*np.ones(M_NZ), np.array([math.log(x) for x in p_l_MAP_NZ]) )
	elif prior == 'Gaussian':
		logL_MAP = np.dot( N_i, np.array([math.log(x) for x in f_i_MAP]) ) + np.dot( M_l_NZ, np.array([math.log(x) for x in p_l_MAP_NZ]) ) - np.dot( alpha*np.ones(M_NZ), np.array([math.pow(x,genGaussPriorPower) for x in p_l_MAP_NZ]) )
	else:
		print("\nERROR - Input prior = '%s' unrecognized, must be one of {'none','Dirichlet','Gaussian'}" % (prior))
		sys.exit(-1)
	fout.write("%15.5e\n" % (logL_MAP))


# Metropolis-Hastings sampling of Bayes posterior distribution to determine uncertainty in the MAP estimate of unbiased probability distribution 

# - masking out all zero-count bins not visited by any simulation for computational efficiency; 
#   accordingly perturbations considered only within bins visited by biased simulations 
#   -> optimal probability estimate for unvisited bins with M_l=0 is p_l=0, so can eliminate these bins from consideration 
#   -> Ref: C. Bartels "Analyzing biased Monte Carlo and molecular dynamics simulations" Chemical Physics Letters 331 446-454 (2000)

# - initializing rng
np.random.seed(seed_MH)

# - performing steps_MH rounds of Metropolis-Hastings sampling
#   -> proposing trial moves as uniform random perturbations on [-maxDpStep_MH, maxDpStep_MH] of randomly selected bin of p_l and renormalizing; 
#      since trial moves are symmetric, this is actually Metropolis sampling which is a specialization of Metropolis-Hastings 
#   -> saving p_l, f_i, and log(L) (to within an additive constant) every saveMod_MH steps as samples from the posterior distribution
print("Metropolis-Hastings sampling of Bayes posterior to determine uncertainty in MAP estimate of probability distribution...\n")

# - initializing sampling matrices
nSamples_MH = int(math.floor(float(steps_MH)/float(saveMod_MH)))
step_MH = np.zeros(nSamples_MH, dtype=np.uint64)
p_l_MH = np.zeros((M,nSamples_MH))
f_i_MH = np.ones((S,nSamples_MH))*float('nan')
logL_MH = np.ones(nSamples_MH)*float('nan')

# - commencing Metropolis-Hastings sampling
maskNZ = (M_l>0)
M_l_NZ = np.array(M_l[maskNZ])
c_il_NZ = np.array(c_il[:,maskNZ])
p_l_NZ = np.array(p_l_MAP[maskNZ])
M_NZ = M_l_NZ.size
f_i = np.array(f_i_MAP)

hitCount=0
for step in range(0,steps_MH):
	
	# -- saving old state
	p_l_NZ_OLD = np.array(p_l_NZ)
	f_i_OLD = np.array(f_i)
	
	# -- proposing a (symmetric) move in {p_l_NZ,f_i}
	idx = np.random.randint(0,M_NZ)

	flag=True
	p_l_NZ_idx_original = p_l_NZ[idx]
	while (flag):
		p_l_NZ[idx] = p_l_NZ[idx] + np.random.uniform(-1.0,1.0)*maxDpStep_MH
		if ( (p_l_NZ[idx] < 0.0) or ( p_l_NZ[idx] > 1.0) ):
			p_l_NZ[idx] = p_l_NZ_idx_original		# proposed move generates invalid p_l, try again
		else:
			flag=False								# proposed move generates valid p_l, proceed
	p_l_NZ = np.divide( p_l_NZ, np.sum(p_l_NZ) )
	f_i = np.divide( 1, np.dot(c_il_NZ,p_l_NZ) )
    
    # -- Metropolis accept/reject
	if prior == 'none':
		teller = math.exp( np.dot( N_i, np.array([math.log(x) for x in np.divide(f_i,f_i_OLD)]) ) + np.dot( M_l_NZ, np.array([math.log(x) for x in np.divide(p_l_NZ,p_l_NZ_OLD)]) ) )
	elif prior == 'Dirichlet':
		teller = math.exp( np.dot( N_i, np.array([math.log(x) for x in np.divide(f_i,f_i_OLD)]) ) + np.dot( ( M_l_NZ + (alpha-1)*np.ones(M_NZ) ), np.array([math.log(x) for x in np.divide(p_l_NZ,p_l_NZ_OLD)]) ) )
	elif prior == 'Gaussian':
		teller = math.exp( np.dot( N_i, np.array([math.log(x) for x in np.divide(f_i,f_i_OLD)]) ) + np.dot( M_l_NZ, np.array([math.log(x) for x in np.divide(p_l_NZ,p_l_NZ_OLD)]) ) - np.dot( alpha*np.ones(M_NZ), (np.array([math.pow(x,genGaussPriorPower) for x in p_l_NZ]) - np.array([math.pow(x,genGaussPriorPower) for x in p_l_NZ_OLD])) ) )
	else:
		print("\nERROR - Input prior = '%s' unrecognized, must be one of {'none','Dirichlet','Gaussian'}" % (prior))
		sys.exit(-1)
	
	if np.random.uniform(0.0,1.0) < teller:		# accept
		hitCount = hitCount + 1
	else:										# reject
		p_l_NZ = np.array(p_l_NZ_OLD)
		f_i = np.array(f_i_OLD)
	
	# -- saving sample to sampling arrays
	if ( (step+1) % saveMod_MH == 0 ):
		jj = (step+1)//saveMod_MH - 1
		step_MH[jj] = step+1
		p_l_MH[maskNZ,jj] = p_l_NZ
		f_i_MH[:,jj] = f_i
		if prior == 'none':
			logL_MH[jj] = np.dot( N_i, np.array([math.log(x) for x in f_i]) ) + np.dot( M_l_NZ, np.array([math.log(x) for x in p_l_NZ]) )
		elif prior == 'Dirichlet':
			logL_MH[jj] = np.dot( N_i, np.array([math.log(x) for x in f_i]) ) + np.dot( M_l_NZ + (alpha-1)*np.ones(M_NZ), np.array([math.log(x) for x in p_l_NZ]) )
		elif prior == 'Gaussian':
			logL_MH[jj] = np.dot( N_i, np.array([math.log(x) for x in f_i]) ) + np.dot( M_l_NZ, np.array([math.log(x) for x in p_l_NZ]) ) - np.dot( alpha*np.ones(M_NZ), np.array([math.pow(x,genGaussPriorPower) for x in p_l_NZ]) )
		else:
			print("\nERROR - Input prior = '%s' unrecognized, must be one of {'none','Dirichlet','Gaussian'}" % (prior))
			sys.exit(-1)
	    
    # -- printing progress to screen
	if ( (step+1) % printMod_MH == 0 ):
		print("\tStep %d of %d, acceptance ratio = %15.5e" % (step+1,steps_MH,float(hitCount)/float(step+1)))
		
print("DONE!\n\n")

# - converting p_l_MH into pdf_l_MH and betaF_l_MH
#   -> mean zeroing betaF; when plotting multiple free energy surfaces this is equivalent to minimizing the mean squared error between the landscapes 
pdf_l_MH = np.zeros(p_l_MH.shape)
betaF_l_MH = np.ones(p_l_MH.shape)*float('nan')
for l in range(0,M):
	sub = ind2sub_RMO(M_k,l)
	binVol = 1
	for d in range(0,dim):
		binVol = binVol*binW[d][sub[d]]
	for jj in range(0,nSamples_MH):
		if p_l_MH[l,jj] > 0:
			pdf_l_MH[l,jj] = p_l_MH[l,jj]/binVol
			betaF_l_MH[l,jj] = -math.log(p_l_MH[l,jj]/binVol)
betaF_l_MH_colMean = np.ones(nSamples_MH)*float('nan')
for jj in range(0,nSamples_MH):
	betaF_l_MH_col = betaF_l_MH[:,jj]
	betaF_l_MH_colMean[jj] = np.mean(betaF_l_MH_col[np.isfinite(betaF_l_MH_col)])
betaF_l_MH = betaF_l_MH - np.matlib.repmat( betaF_l_MH_colMean, betaF_l_MH.shape[0], 1 )

# - writing to file nSamples_MH Metropolis-Hastings samples from Bayes posterior of p_l, pdf_l, and betaF_l over l=1..M bins of dim-dimensional rectilinear histogram grid as row vectors in row major order (last index changes fastest) 
with open("p_MH.txt","w") as fout:
	for jj in range(0,nSamples_MH):
		for l in range(0,M):
			fout.write("%15.5e" % (p_l_MH[l,jj]))
		fout.write("\n")

with open("pdf_MH.txt","w") as fout:
	for jj in range(0,nSamples_MH):
		for l in range(0,M):
			fout.write("%15.5e" % (pdf_l_MH[l,jj]))
		fout.write("\n")

with open("betaF_MH.txt","w") as fout:
	for jj in range(0,nSamples_MH):
		for l in range(0,M):
			fout.write("%15.5e" % (betaF_l_MH[l,jj]))
		fout.write("\n")
		
# - writing to file nSamples_MH Metropolis-Hastings samples from Bayes posterior of f_i = Z/Z_i = ratio of unbiased partition function to that of biased simulation i,for i=1..S biased simulations 
with open("f_MH.txt","w") as fout:
	for jj in range(0,nSamples_MH):
		for i in range(0,S):
			fout.write("%15.5e" % (f_i_MH[i,jj]))
		fout.write("\n")

# - writing to file log of Bayes posterior up to an additive constant for the nSamples_MH Metropolis-Hastings samples from the posterior 
with open("logL_MH.txt","w") as fout:
	for jj in range(0,nSamples_MH):
		fout.write("%15.5e\n" % (logL_MH[jj]))

# - writing to file associated MH step of the nSamples_MH Metropolis-Hastings samples from Bayes posterior  
with open("step_MH.txt","w") as fout:
	for jj in range(0,nSamples_MH):
		fout.write("%15d\n" % (step_MH[jj]))





#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# OVERVIEW
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Purpose: the purpose of this script is to conduct LRT within DFE type and estimate parameter uncertainties using the Godambe methods for composite likelihoods
# Output:
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Run as: python godambe.py 
#
# ~~~Required base arguments~~~
# --prefix [prefix to use for output/population(s) name(s)] 
# --fs [fs file for selected sites]
# --dfe_model [argument for DFE model]
# --cache1d [name of pickled cached spectra for 1D case] 
# --cache2d [name of pickled cached spectra for 2D case] 
# --mut [mutation rate] 
# --L_neut [effective sequence length for neutral compartment]
# --L_sel [effective sequence length for selected compartment] 
# --theta_neut [inferred theta for neutral demographic model]
#
# ~~~Needed to make bootstrap replicates~~~
# --make_boot [boolean that determines whether to make the bootstrapped datasets for the given FS file]
# --vcf [required if using the --make_boot option]
# --popfile [required if using the --make_boot option]
#
# ~~~Needed to perform LRT~~~
# --lrt [boolean that determines whether to conduct LRT]
# --complex_model [argument to specify whether complex model is the addition of pneu or ppos]
# --simple_popt [best-fit parameter estimates for the simpler model]
# --complex_popt [best-fit parameter estimates for the more complex model]
#
# ~~~Needed to estimate uncertainties~~~
# --uncert [boolean that determines whether to estimate model uncertainties]
# --best_fit_model [argument to specify whether best-fit model is base, pneu, or ppos variations of --dfe_model]
# --best_fit_popt [specify the parameters of the best-fit model]
#
# ~~~Needed to compute CL-AIC metric~~~
# --claic [boolean that determines whether to compute the CL-AIC]
# --best_fit_model [argument to specify whether best-fit model is base, pneu, or ppos variations of --dfe_model]
# --best_fit_popt [specify the parameters of the best-fit model]
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#--------------------------
# Import required packages
#--------------------------
import dadi
import random
import numpy as np
import matplotlib.pylab as pylab
import matplotlib.pyplot as pyplot
import argparse
import dadi.DFE as DFE
import dill as pickle
import glob
import scipy.stats.distributions as ssd

#-------------------------------------------------
# Loading the data and reading command line input
#-------------------------------------------------

def main():

	print('---------------------------------------------------------------------------')
	print('INITIALIZING')
	print('---------------------------------------------------------------------------')

	print('Parsing command-line arguments...')

	# Parse args
	args = parse_args()

	# Required arguments

	# Create variable to store the population name
	prefix = args.prefix
	# Read in the input frequency spectrum
	fs = dadi.Spectrum.from_file(args.fs)
	# Determine the DFE model type to fit
	dfe_model = args.dfe_model
	# Load the 1D and 2D cached spectra
	cache1d = pickle.load(open(args.cache1d,'rb'))
	cache2d = pickle.load(open(args.cache2d,'rb'))
	# Placeholder to store bad gamma values
	bad_gamma_union = []

	print('Checking for NaN values in the 1D cached spectra...')
	# Check if the 1D cache contains NaN values
	if np.isnan(cache1d.spectra).any():
		print('Found NaN values...')
		# Determine which of the 1D cached 2D SFS (which are ordered by gamma) contain the NaN values
		bad_gamma_1d = cache1d.gammas[np.isnan(cache1d.spectra).any(axis=(1,2))]
		print('Offending gamma value(s): ' + str(bad_gamma_1d))
		# Add these bad gammas to the list
		bad_gamma_union = bad_gamma_union + bad_gamma_1d.tolist()
	else:
		print('None found!')

	print('Checking for NaN values in the 2D cached spectra...')
	# Check if the 2D cache contains NaN values
	if np.isnan(cache2d.spectra).any():
		print('Found more NaN values...')
		# Determine which of the 2D cached 2D SFS (which are ordered by a grid of gamma) contain the NaN values along the first dimension (0)
		bad_gamma_2d_dim1 = cache2d.gammas[np.isnan(cache2d.spectra).any(axis=(1,2,3))]
		print('Offending gamma value(s) in dimension 1: ' + str(bad_gamma_2d_dim1))
		# Add these bad gammas to the list
		bad_gamma_union = bad_gamma_union + bad_gamma_2d_dim1.tolist()
	else:
		print('None found!')

	# Find the unique union of bad gammas across 1D and 2D caches
	bad_gamma_union = np.unique(bad_gamma_union)
	
	if len(bad_gamma_union > 0):

		# Set the correction value to be the "maximum" of these negative values
		correction = max(bad_gamma_union)

		print('Correcting cached frequency spectra to exclude negative gammas equal to or exceeding ' + str(correction))

		print('Correcting 1D cache...')

		# Find the index of the cutoff
		cutoff_idx = min(np.where(cache1d.gammas >= correction + 1)[0].tolist())

		# Only keep indices that are less than (in abs value) the cutoff
		cache1d.spectra = cache1d.spectra[cutoff_idx:]

		# Perform the same operation to the cache.gammas and cache.neg_gammas attributes
		cache1d.gammas = cache1d.gammas[cutoff_idx:]
		cache1d.neg_gammas = cache1d.neg_gammas[cutoff_idx:]
		
		print('Correcting 2D cache...')

		# Find the index of the cutoff
		cutoff_idx = min(np.where(cache2d.gammas >= correction + 1)[0].tolist())

		# Only keep indices that are less than (in abs value) the cutoff
		cache2d.spectra = cache2d.spectra[cutoff_idx:, cutoff_idx:]

		# Perform the same operation to the cache.gammas and cache.neg_gammas attributes
		cache2d.gammas = cache2d.gammas[cutoff_idx:]
		cache2d.neg_gammas = cache2d.neg_gammas[cutoff_idx:]

	# Create a variable to hold the specified mutation rate
	mu = float(args.mut)
	# Create a variable to hold the specified effective sequence length for neutral region
	L_neut = int(args.L_neut)
	# Create a variable to hold the specified effective sequence length for selected region
	L_sel = int(args.L_sel)
	# Parameter to hold inferred theta for the neutral case
	theta_neut = int(args.theta_neut)

	# Arguments for making bootstrap replicates

	# Store a boolean that determines whether we should make bootstrapped datasets for the given data
	make_boot = args.make_boot
	# Store the required VCF for the --make_boot option
	vcf = args.vcf
	# Store the popfile required for the --make_boot option
	popfile = args.popfile

	# Arguments for conducting likelihood ratio tests

	# Store a boolean that determines whether to conduct LRT
	lrt = args.lrt
	# Create a variable that stores the type of the more complex model (either pneu or ppos)
	complex_model = args.complex_model
	# Store a list of the simple model (specified above with --dfe_model) parameter estimates
	simple_popt = [float(p) for p in args.simple_popt]
	# Store a list of the complex model (i.e., dfe_model + complex_model) parameter estimates
	complex_popt = [float(p) for p in args.complex_popt]

	# Arguments for estimating parameter uncertainties for a given model

	# Store a boolean that determines whether to estimate parameter uncertainties
	uncert = args.uncert
	# Create a variable that specifies the complexity of the best-fit DFE model (as either base, pneu, or ppos)
	best_fit_model = args.best_fit_model
	# Store a list of the parameter estimates for the best-fit model
	best_fit_popt = [float(p) for p in args.best_fit_popt]

	# Additional arguments for computing the CL-AIC (also used best_fit_model and best_fit_popt from above)

	# Store a boolean that determines whether to compute CL-AIC
	claic = args.claic

	#------------------------
	# Initializing variables
	#------------------------

	# Need to determine what the sample size is based on the input fs
	ns = fs.sample_sizes
	print("Sample sizes:")
	print(ns)

	# Need to define the grid points for optimization
	pts = [x for x in range(0, len(ns)*10, 10)] + (max(ns)*2)
	print("Grid points:")
	print(pts)

	# Need to determine pop_ids based on input fs
	pop_ids = fs.pop_ids
	print("Population IDs:")
	print(pop_ids)

	# Need to compute theta for the selected sites (I had previously been computing this incorrectly)
	theta_sel = theta_neut * (L_sel/L_neut)
	print("Theta for selected sites:")
	print(theta_sel)

	#---------------------------------
	# Making (or loading) bootstraps
	#---------------------------------

	print('---------------------------------------------------------------------------')
	print('BOOTSTRAPPING')
	print('---------------------------------------------------------------------------')

	# Create a directory to store the bootstrapped FS
	boots_dir = 'bootstrap_replicates/'

	# Make the boostraps if specified
	if make_boot:
		print('Making bootstrap FS...')

		# Read in the input vcf and construct dadi dictionary object
		dd = dadi.Misc.make_data_dict_vcf(vcf, popfile)
		# Divide the genome into 2 Mbp sized chunks
		chunks = dadi.Misc.fragment_data_dict(dd, chunk_size=2e6)
		# Create 100 bootstrap replicates by sampling with replacement from 2 Mbp chunks
		boots = dadi.Misc.bootstraps_from_dd_chunks(chunks, Nboot=100, pop_ids=pop_ids, projections=ns, polarized=True)

		# Index bootstrap replicates
		index = 1
		for boot in boots:
			print('Writing bootstrap replicate ' + str(index) + '...')
			# Write each boostrap to an indexed file
			boot_name = boots_dir + str(index) + '_boot_' + prefix + '.fs'
			boot.to_file(boot_name)
			# Increase the index
			index = index + 1

		print ('Done.')

	print('Reading in bootstrap FS...')

	# Load the saved bootstrap replicates for the given dataset
	boots_glob = boots_dir + "*" + '_boot_' + prefix + '.fs'
	boots_fids = glob.glob(boots_glob)
	boots_fs = [dadi.Spectrum.from_file(fid) for fid in boots_fids]

	# Generate a sequence of relative values of theta (compared to original data) to use for Godambe functions
	boot_theta_adjusts = [boot.sum()/fs.sum() for boot in boots_fs]

	#---------------------------------------------------------------------------------------------------------------
	# Conducting LRT within one of the three DFE types (biv_lognormal_shared, biv_lognormal_ind, lognormal_mixture)
	#---------------------------------------------------------------------------------------------------------------

	if lrt:

		print('---------------------------------------------------------------------------')
		print('LIKELIHOOD RATIO TEST')
		print('---------------------------------------------------------------------------')

		print('Defining simple and complex DFE functions...')

		# Determine which base DFE model to use, followed by which complex DFE model to use
		if dfe_model=='biv_lognormal_shared':
			# Store the simple model function
			simple_func, simple_dist = biv_lognormal_shared(cache2d)
			# Store the complex model function
			if complex_model=='pneu':
				complex_func, complex_dist, nested_indices, fixed_indices = biv_lognormal_shared_pneu(cache2d)
			elif complex_model=='ppos':
				complex_func, complex_dist, nested_indices, fixed_indices = biv_lognormal_shared_ppos(cache2d)
		if dfe_model=='biv_lognormal_ind':
			# Store the simple model function
			simple_func, simple_dist = biv_lognormal_ind(cache2d)
			# Store the complex model function
			if complex_model=='pneu':
				complex_func, complex_dist, nested_indices, fixed_indices = biv_lognormal_ind_pneu(cache2d)
			elif complex_model=='ppos':
				complex_func, complex_dist, nested_indices, fixed_indices = biv_lognormal_ind_ppos(cache2d)
		if dfe_model=='mixture_lognormal':
			# Store the simple model function
			simple_func, simple_dist, fixed_indices = mixture_lognormal()
			# Store the complex model function
			if complex_model=='pneu':
				complex_func, complex_dist, nested_indices, fixed_indices = mixture_lognormal_pneu()
			elif complex_model=='ppos':
				complex_func, complex_dist, nested_indices, fixed_indices = mixture_lognormal_ppos()

		print('Computing likelihood for simple model with simple parameters...')

		# Use the fixed_indices to get the values of the fixed parameters in the more complex model
		fixed_vals = [complex_popt[idx] for idx in fixed_indices]

		# Generate expected FS under simple model
		if len(simple_dist) > 1:
			# If more than one simple_dist is specified, then provide arguments to the cache.mixture function
			s_model = simple_func(simple_popt, ns, cache1d, cache2d, simple_dist[0], simple_dist[1], theta_sel, pts)
		elif len(simple_dist) == 1:
			# If only one simple_dist is specified, then provide arguments to the cache.integrate function
			s_model = simple_func(simple_popt, ns, simple_dist[0], theta_sel, pts)

		# Compute likelihood under simple model
		s_ll = dadi.Inference.ll_multinom(s_model, fs)
		print(s_ll)

		print('Computing likelihood for complex model with complex parameters...')

		# Generate expected FS under the complex model
		if len(complex_dist) > 1:
			# If more than one complex_dist is specified, then provide arguments to the cache.mixture function
			c_model = complex_func(complex_popt, ns, cache1d, cache2d, complex_dist[0], complex_dist[1], theta_sel, pts)
		elif len(complex_dist) == 1:
			# If only one complex_dist is specified, then provide arguments to the cache.integrate function
			c_model = complex_func(complex_popt, ns, complex_dist[0], theta_sel, pts)

		# Compute likelihood under complex model
		c_ll = dadi.Inference.ll_multinom(c_model, fs)
		print(c_ll)

		print('Creating parameters for complex model using best-fit estimates from simple model...')

		# Define new parameter list (inputting simple model estimates into complex model parameters and setting nested indices to 0)
		p_lrt = simple_popt
		for i in nested_indices:
			p_lrt.insert(i, 0)

		# Delete any fixed_indices, which will be added back (using the fixed values) in the new_complex_func function
		p_lrt = np.delete(p_lrt, fixed_indices)
		
		print('These are the simple model estimates input into the complex model parameter list, with fixed parameters removed:')
		print(p_lrt)

		# The nested indices used in the call to Godambe.LRT_adjust should exclude any fixed indices
		new_nested_indices = [elem for elem in nested_indices if elem not in fixed_indices]
		
		print('These are the nested model indices:')
		print(new_nested_indices)

		print('Computing Godambe adjustment for LRT...')

		# Create a new dadi function for the complex model that includes the extra arguments needed for DFE analysis
		def new_complex_func(params, ns, pts):
			# First, insert placeholder elements (i.e., zeros) for parameters that are to be fixed
			for i in fixed_indices:
				params = np.insert(params, i, 0)
			# Then, populate these placeholders with the corresponding fixed values
			for i in range(0, len(fixed_indices)):
				params[fixed_indices[i]] = fixed_vals[i]
			print('These are the simple model estimates input into the complex model parameter list, with fixed parameters included and nested indices set to 0:')
			print(params)
			# Generate expected FS under complex model
			if len(complex_dist) > 1:
				# If more than one complex_dist is specified, then provide arguments to the cache.mixture function
				return complex_func(params, ns, cache1d, cache2d, complex_dist[0], complex_dist[1], theta_sel, pts)
			elif len(complex_dist) == 1:
				# If only one complex_dist is specified, then provide arguments to the cache.integrate function
				return complex_func(params, ns, complex_dist[0], theta_sel, pts)

		# Loop through different possible step sizes to ensure uncertainty estimates are stable
		for eps in [0.01, 0.001, 0.0001]:
		
			# Compute the adjustment using p_lrt and the bootstraps created above
			adj = dadi.Godambe.LRT_adjust(new_complex_func, pts, boots_fs, p_lrt, fs, nested_indices=new_nested_indices, multinom=False, boot_theta_adjusts=boot_theta_adjusts, eps=eps)
			# Compute the adjusted LRT test statistic
			D_adj = adj*2*(c_ll - s_ll)
			# Compute the p-value (these weights assumes the models differ by only a single free parameter)
			pval = dadi.Godambe.sum_chi2_ppf(D_adj, weights=(0.5, 0.5))

			print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
			print("Using step size = " + str(eps))
			print ('Godambe adjustment value:')
			print(adj)
			print('Adjusted LRT test statistic:')
			print(D_adj)
			print('P-value for rejecting simpler model:')
			print(pval)
			if pval < 0.05:
				print('The complex model provides a significantly better fit')
			else:
				print('The simpler model provides a significantly better fit')
			print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

	#---------------------------------------------------------------------------------------------------------------
	# Estimating parameter uncertainties for the given model
	#---------------------------------------------------------------------------------------------------------------

	if uncert:

		print('---------------------------------------------------------------------------')
		print('PARAMETER UNCERTAINTIES')
		print('---------------------------------------------------------------------------')

		print('Defining best-fit DFE function...')

		# Variable to initialize an empty array for the fixed_indices (in the event that we're estimating parameters for a base model)
		fixed_indices=[]

		# Determine which base DFE model to use, followed by which complex DFE model to use
		if dfe_model=='biv_lognormal_shared':
			# Determine whether we want the base model, the model with neutral point mass, or the model with positive point mass
			if best_fit_model=='base':
				dfe_func, dfe_dist = biv_lognormal_shared(cache2d)
			elif best_fit_model=='pneu':
				dfe_func, dfe_dist, nested_indices, fixed_indices = biv_lognormal_shared_pneu(cache2d)
			elif best_fit_model=='ppos':
				dfe_func, dfe_dist, nested_indices, fixed_indices = biv_lognormal_shared_ppos(cache2d)
		if dfe_model=='biv_lognormal_ind':
			# Determine whether we want the base model, the model with neutral point mass, or the model with positive point mass
			if best_fit_model=='base':
				dfe_func, dfe_dist = biv_lognormal_ind(cache2d)
			elif best_fit_model=='pneu':
				dfe_func, dfe_dist, nested_indices, fixed_indices = biv_lognormal_ind_pneu(cache2d)
			elif best_fit_model=='ppos':
				dfe_func, dfe_dist, nested_indices, fixed_indices = biv_lognormal_ind_ppos(cache2d)
		if dfe_model=='mixture_lognormal':
			# Determine whether we want the base model, the model with neutral point mass, or the model with positive point mass
			if best_fit_model=='base':
				dfe_func, dfe_dist, fixed_indices = mixture_lognormal()
			elif best_fit_model=='pneu':
				dfe_func, dfe_dist, nested_indices, fixed_indices = mixture_lognormal_pneu()
			elif best_fit_model=='ppos':
				dfe_func, dfe_dist, nested_indices, fixed_indices = mixture_lognormal_ppos()

		# Best-fit parameter estimates
		popt = best_fit_popt 
		# Use the fixed_indices to get the values of the fixed parameters 
		fixed_vals = [popt[idx] for idx in fixed_indices]

		print('Best-fit parameter estimates (including p_misid):')
		print(popt)
		print('Indices of fixed parameters:')
		print(fixed_indices)
		print('Values of fixed parameters:')
		print(fixed_vals)

		# Exclude the fixed parameter estimates from the list passed to the Godambe function
		popt = np.delete(popt, fixed_indices)
		print('These are the best-fit parameter estimates, with fixed parameters removed:')
		print(popt)

		# Create a new dadi function for the DFE model that includes the extra arguments needed for DFE analysis and allows us to fix certain parameters
		def new_dfe_func(params, ns, pts):
			# First, insert placeholder elements (i.e., zeros) for parameters that are to be fixed
			for i in fixed_indices:
				params = np.insert(params, i, 0)
			# Then, populate these placeholders with the corresponding fixed values
			for i in range(0, len(fixed_indices)):
				params[fixed_indices[i]] = fixed_vals[i]
			# Generate expected FS under the best-fit model
			if len(dfe_dist) > 1:
				# If more than one dfe_dist is specified, then provide arguments to the cache.mixture function
				return dfe_func(params, ns, cache1d, cache2d, dfe_dist[0], dfe_dist[1], theta_sel, pts)
			elif len(dfe_dist) == 1:
				# If only one dfe_dist is specified, then provide arguments to the cache.integrate function
				return dfe_func(params, ns, dfe_dist[0], theta_sel, pts)
		
		print("Estimating parameter uncertainties...")
		# Loop through different possible step sizes to ensure uncertainty estimates are stable
		for eps in [0.01, 0.001, 0.0001]:

			# Function returns standard deviation of parameter values (including theta– listed last) along with the full GIM to use in propogating uncertainties
			uncert_est, GIM, H = dadi.Godambe.GIM_uncert(new_dfe_func, pts, boots_fs, popt, fs, multinom=False, eps=eps, return_GIM=True, boot_theta_adjusts=boot_theta_adjusts)

			print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
			print("Using step size = " + str(eps))
			print("Standard deviations of *free* parameter values:")
			print(uncert_est)
			print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

	#---------------------------------------------------------------------------------------------------------------
	# Computing the composite likelihood of AIC for the given model
	#---------------------------------------------------------------------------------------------------------------

	if claic:

		print('---------------------------------------------------------------------------')
		print('COMPUTING CL-AIC')
		print('---------------------------------------------------------------------------')

		print('Defining best-fit DFE function...')

		# Variable to initialize an empty array for the fixed_indices (in the event that we're estimating parameters for a base model)
		fixed_indices=[]

		# Determine which base DFE model to use, followed by which complex DFE model to use
		if dfe_model=='biv_lognormal_shared':
			# Determine whether we want the base model, the model with neutral point mass, or the model with positive point mass
			if best_fit_model=='base':
				dfe_func, dfe_dist = biv_lognormal_shared(cache2d)
			elif best_fit_model=='pneu':
				dfe_func, dfe_dist, nested_indices, fixed_indices = biv_lognormal_shared_pneu(cache2d)
			elif best_fit_model=='ppos':
				dfe_func, dfe_dist, nested_indices, fixed_indices = biv_lognormal_shared_ppos(cache2d)
		if dfe_model=='biv_lognormal_ind':
			# Determine whether we want the base model, the model with neutral point mass, or the model with positive point mass
			if best_fit_model=='base':
				dfe_func, dfe_dist = biv_lognormal_ind(cache2d)
			elif best_fit_model=='pneu':
				dfe_func, dfe_dist, nested_indices, fixed_indices = biv_lognormal_ind_pneu(cache2d)
			elif best_fit_model=='ppos':
				dfe_func, dfe_dist, nested_indices, fixed_indices = biv_lognormal_ind_ppos(cache2d)
		if dfe_model=='mixture_lognormal':
			# Determine whether we want the base model, the model with neutral point mass, or the model with positive point mass
			if best_fit_model=='base':
				dfe_func, dfe_dist, fixed_indices = mixture_lognormal()
			elif best_fit_model=='pneu':
				dfe_func, dfe_dist, nested_indices, fixed_indices = mixture_lognormal_pneu()
			elif best_fit_model=='ppos':
				dfe_func, dfe_dist, nested_indices, fixed_indices = mixture_lognormal_ppos()

		# Best-fit parameter estimates
		popt = best_fit_popt
		# Use the fixed_indices to get the values of the fixed parameters 
		fixed_vals = [popt[idx] for idx in fixed_indices]

		print('Best-fit parameter estimates (including p_misid):')
		print(popt)
		print('Indices of fixed parameters:')
		print(fixed_indices)
		print('Values of fixed parameters:')
		print(fixed_vals)

		print('Computing likelihood for the best-fit model with the best-fit parameter estimates...')

		# Generate expected FS under given model
		if len(dfe_dist) > 1:
			# If more than one dfe_dist is specified, then provide arguments to the cache.mixture function
			model = dfe_func(popt, ns, cache1d, cache2d, dfe_dist[0], dfe_dist[1], theta_sel, pts)
		elif len(dfe_dist) == 1:
			# If only one dfe_dist is specified, then provide arguments to the cache.integrate function
			model = dfe_func(popt, ns, dfe_dist[0], theta_sel, pts)

		# Compute likelihood under this model
		ll = dadi.Inference.ll_multinom(model, fs)
		print(ll)

		# Exclude the fixed parameter estimates from the list passed to the Godambe function
		popt = np.delete(popt, fixed_indices)
		print('These are the best-fit parameter estimates, with fixed parameters removed:')
		print(popt)

		# Create a new dadi function for the DFE model that includes the extra arguments needed for DFE analysis and allows us to fix certain parameters
		def new_dfe_func(params, ns, pts):
			# First, insert placeholder elements (i.e., zeros) for parameters that are to be fixed
			for i in fixed_indices:
				params = np.insert(params, i, 0)
			# Then, populate these placeholders with the corresponding fixed values
			for i in range(0, len(fixed_indices)):
				params[fixed_indices[i]] = fixed_vals[i]
			# Generate expected FS under the best-fit model
			if len(dfe_dist) > 1:
				# If more than one dfe_dist is specified, then provide arguments to the cache.mixture function
				return dfe_func(params, ns, cache1d, cache2d, dfe_dist[0], dfe_dist[1], theta_sel, pts)
			elif len(dfe_dist) == 1:
				# If only one dfe_dist is specified, then provide arguments to the cache.integrate function
				return dfe_func(params, ns, dfe_dist[0], theta_sel, pts)
		
		print("Obtaining Godambe information and Hessian matrices...")
		# Loop through different possible step sizes to ensure Godambe results are stable
		for eps in [0.01, 0.001, 0.0001]:
		#for eps in [0.001, 0.0001]:
			# Function returns the GIM, H matrix, J matrix and cU column vector to use in calculating the CL-AIC statistic
			GIM, H, J, cU = dadi.Godambe.get_godambe(new_dfe_func, pts, boots_fs, popt, fs, eps=eps, boot_theta_adjusts=boot_theta_adjusts)
			print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
			print("Using step size = " + str(eps))
			print('Computing CL-AIC as [2 * np.trace(np.dot(J, np.linalg.inv(H))) - 2 * ll]:')
			CLAIC = 2 * np.trace(np.dot(J, np.linalg.inv(H))) - 2 * ll
			print(CLAIC)
			print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
			print("Note: The smaller the AIC or CLAIC score is, the better the model fits the observed data")

#----------------------------------
# Defining functions
#----------------------------------
def biv_lognormal_shared(cache2d):
	"""
	Prepares all necessary variables to fit a bivariate lognormal distribution (with shared mu and sigma) to the DFE
	params = mu, sigma, rho
	"""
	print("Initializing bivariate lognormal joint DFE with shared parameters...")

	# Define the DFE function used (either cache.integrate for bivariate prob. dist or cache.mixture for mixture models)
	dfe_func = cache2d.integrate

	print("Adding a parameter for ancestral state misidentification")
	# Add a parameter for ancestral state misidentification
	dfe_func = dadi.Numerics.make_anc_state_misid_func(dfe_func)

	print("Defining selection distribution")
	# Define which selection distribution to use (a single 2D PDF for the bivariate prob. dist or a 1D and 2D PDF for the mixture model)
	dfe_dist = DFE.PDFs.biv_lognormal

	# Return relevant arguments (note that dfe_dist is only returned in this way to allow automatic calculation of the model FS in the main function)
	return dfe_func, [dfe_dist]

def biv_lognormal_shared_pneu(cache2d):
	"""
	Prepares all necessary variables to fit a bivariate lognormal distribution (with shared mu and sigma) plus a symmetric point mass of neutrality to the DFE
	params = mu, sigma, rho, pneu, gammaneu
	* pneu defines the neutral proportion
	* gammaneu defines the neutral selection coefficient (fixed to 0)
	"""
	print("Initializing bivariate lognormal joint DFE with shared parameters *and* a symmetric point mass of neutrality...")

	# Define the DFE function to use (here we use the integration method for a symmetric point mass of positive selection to model a neutral component)
	dfe_func = cache2d.integrate_symmetric_point_pos

	print("Adding a parameter for ancestral state misidentification")
	# Add a parameter for ancestral state misidentification
	dfe_func = dadi.Numerics.make_anc_state_misid_func(dfe_func)

	print("Defining selection distribution")
	# Define which selection distribution to use (a single 2D PDF for the bivariate prob. dist or a 1D and 2D PDF for the mixture model)
	dfe_dist = DFE.PDFs.biv_lognormal

	# Define the nested and fixed indices (starting at 0)
	nested_indices = [3, 4]
	fixed_indices = [4]

	# Return relevant arguments (note that dfe_dist is only returned in this way to allow automatic calculation of the model FS in the main function)
	return dfe_func, [dfe_dist], nested_indices, fixed_indices

def biv_lognormal_shared_ppos(cache2d):
	"""
	Prepares all necessary variables to fit a bivariate lognormal distribution (with shared mu and sigma) plus a symmetric point mass of positive selection to the DFE
	params = mu, sigma, rho, ppos, gammapos
	* ppos defines the positive proportion
	* gammapos defines the positive selection coefficient (fixed at the value specified on the command line)
	"""
	print("Initializing bivariate lognormal joint DFE with shared parameters *and* a symmetric point mass of positive selection...")

	# Define the DFE function to use (here we use the integration method for a symmetric point mass of positive selection)
	dfe_func = cache2d.integrate_symmetric_point_pos

	print("Adding a parameter for ancestral state misidentification")
	# Add a parameter for ancestral state misidentification
	dfe_func = dadi.Numerics.make_anc_state_misid_func(dfe_func)

	print("Defining selection distribution")
	# Define which selection distribution to use (a single 2D PDF for the bivariate prob. dist or a 1D and 2D PDF for the mixture model)
	dfe_dist = DFE.PDFs.biv_lognormal

	# Define the nested and fixed indices (starting at 0)
	nested_indices = [3, 4]
	fixed_indices = [4]

	# Return relevant arguments (note that dfe_dist is only returned in this way to allow automatic calculation of the model FS in the main function)
	return dfe_func, [dfe_dist], nested_indices, fixed_indices

def biv_lognormal_ind(cache2d):
	"""
	Prepares all necessary variables to fit a bivariate lognormal distribution (with distinct mu and sigma) to the DFE
	params = mu1, mu2, sigma1, sigma2, rho
	"""
	print("Initializing bivariate lognormal joint DFE with distinct parameters...")

	# Define the DFE function used (either cache.integrate for bivariate prob. dist or cache.mixture for mixture models)
	dfe_func = cache2d.integrate

	print("Adding a parameter for ancestral state misidentification")
	# Add a parameter for ancestral state misidentification
	dfe_func = dadi.Numerics.make_anc_state_misid_func(dfe_func)

	print("Defining selection distribution")
	# Define which selection distribution to use (a single 2D PDF for the bivariate prob. dist or a 1D and 2D PDF for the mixture model)
	dfe_dist = DFE.PDFs.biv_lognormal

	# Return relevant arguments (note that dfe_dist is only returned in this way to allow automatic calculation of the model FS in the main function)
	return dfe_func, [dfe_dist]

def biv_lognormal_ind_pneu(cache2d):
	"""
	Prepares all necessary variables to fit a bivariate lognormal distribution (with distinct mu and sigma) and a symmetric point mass of neutrality to the DFE
	params = mu1, mu2, sigma1, sigma2, rho, pneu, gammaneu
	* pneu defines the neutral proportion
	* gammaneu defines the neutral selection coefficient (fixed to 0)
	"""
	print("Initializing bivariate lognormal joint DFE with distinct parameters *and* a symmetric point mass of neutrality...")

	# Define the DFE function to use (here we use the integration method for a symmetric point mass of positive selection to model a neutral component)
	dfe_func = cache2d.integrate_symmetric_point_pos

	print("Adding a parameter for ancestral state misidentification")
	# Add a parameter for ancestral state misidentification
	dfe_func = dadi.Numerics.make_anc_state_misid_func(dfe_func)

	print("Defining selection distribution")
	# Define which selection distribution to use (a single 2D PDF for the bivariate prob. dist or a 1D and 2D PDF for the mixture model)
	dfe_dist = DFE.PDFs.biv_lognormal

	# Define the nested and fixed indices (starting at 0)
	nested_indices = [5, 6]
	fixed_indices = [6]

	# Return relevant arguments (note that dfe_dist is only returned in this way to allow automatic calculation of the model FS in the main function)
	return dfe_func, [dfe_dist], nested_indices, fixed_indices

def biv_lognormal_ind_ppos(cache2d):
	"""
	Prepares all necessary variables to fit a bivariate lognormal distribution (with distinct mu and sigma) and a symmetric point mass of positive selection to the DFE
	params = mu1, mu2, sigma1, sigma2, rho, ppos, gammapos
	* ppos defines the neutral proportion
	* gammapos defines the neutral selection coefficient (specified on the command line)
	"""
	print("Initializing bivariate lognormal joint DFE with distinct parameters *and* a symmetric point mass of positive selection...")

	# Define the DFE function to use (here we use the integration method for a symmetric point mass of positive selection)
	dfe_func = cache2d.integrate_symmetric_point_pos

	print("Adding a parameter for ancestral state misidentification")
	# Add a parameter for ancestral state misidentification
	dfe_func = dadi.Numerics.make_anc_state_misid_func(dfe_func)

	print("Defining selection distribution")
	# Define which selection distribution to use (a single 2D PDF for the bivariate prob. dist or a 1D and 2D PDF for the mixture model)
	dfe_dist = DFE.PDFs.biv_lognormal

	# Define the nested and fixed indices (starting at 0)
	nested_indices = [5, 6]
	fixed_indices = [6]

	# Return relevant arguments (note that dfe_dist is only returned in this way to allow automatic calculation of the model FS in the main function)
	return dfe_func, [dfe_dist], nested_indices, fixed_indices

def mixture_lognormal():
	"""
	Prepares all necessary variables to fit a mixture model to the DFE
	> This is a mixture of a 1D distribution (representing perfectly correlated fitness effects) and a 2D distribution (representing uncorrelated fitness effects)
	> The 1D distribution is a lognormal distribution (mu, sigma)
	> The 2D distribution is a bivariate lognormal distribution with shared parameters and rho = 0 (mu, sigma, rho=0)
	> The lognormal parameters mu and sigma are shared between the 1D and 2D component
	> The w parameter is the weight of the 2D distribution in the mixture model (interpreted as a DFE correlation)
	params = mu, sigma, rho, w
	"""
	print("Initializing lognormal mixture model...")

	# Define the DFE function used (either cache.integrate for bivariate prob. dist or cache.mixture for mixture models)
	dfe_func = DFE.mixture

	print("Adding a parameter for ancestral state misidentification")
	# Add a parameter for ancestral state misidentification
	dfe_func = dadi.Numerics.make_anc_state_misid_func(dfe_func)

	print("Defining selection distribution")
	# Define which selection distribution to use (a single 2D PDF for the bivariate prob. dist or a 1D and 2D PDF for the mixture model)
	dfe_dist_1d = DFE.PDFs.lognormal
	dfe_dist_2d = DFE.PDFs.biv_lognormal
	# Define the DFE dist list for the purposes of post-optimization model calculation
	dfe_dist = [dfe_dist_1d, dfe_dist_2d]

	# Define the fixed indices (starting at 0)
	fixed_indices = [2]

	# Return relevant arguments (note that dfe_dist is only returned in this way to allow automatic calculation of the model FS in the main function)
	return dfe_func, dfe_dist, fixed_indices

def mixture_lognormal_pneu():
	"""
	Prepares all necessary variables to fit a mixture model to the DFE with a symmetric point mass of neutrality
	> This is a mixture of a 1D distribution (representing perfectly correlated fitness effects) and a 2D distribution (representing uncorrelated fitness effects)
	> The 1D distribution is a lognormal distribution (mu, sigma)
	> The 2D distribution is a bivariate lognormal distribution with shared parameters and rho = 0 (mu, sigma, rho=0)
	> The lognormal parameters mu and sigma are shared between the 1D and 2D component
	> The w parameter is the weight of the 2D distribution in the mixture model (interpreted as a DFE correlation)
	params = mu, sigma, rho, pneu, gammaneu, w
	"""
	print("Initializing lognormal mixture model with a symmetric point mass of neutrality...")

	# Define the DFE function used (either cache.integrate for bivariate prob. dist or cache.mixture for mixture models)
	dfe_func = DFE.mixture_symmetric_point_pos # note that the parameter order is slightly different for this case

	print("Adding a parameter for ancestral state misidentification")
	# Add a parameter for ancestral state misidentification
	dfe_func = dadi.Numerics.make_anc_state_misid_func(dfe_func)

	print("Defining selection distribution")
	# Define which selection distribution to use (a single 2D PDF for the bivariate prob. dist or a 1D and 2D PDF for the mixture model)
	dfe_dist_1d = DFE.PDFs.lognormal
	dfe_dist_2d = DFE.PDFs.biv_lognormal
	# Define the DFE dist list for the purposes of post-optimization model calculation
	dfe_dist = [dfe_dist_1d, dfe_dist_2d]

	# Define the nested and fixed indices (starting at 0)
	nested_indices = [3, 4]
	fixed_indices = [2, 4]

	# Return relevant arguments (note that dfe_dist is only returned in this way to allow automatic calculation of the model FS in the main function)
	return dfe_func, dfe_dist, nested_indices, fixed_indices

def mixture_lognormal_ppos():
	"""
	Prepares all necessary variables to fit a mixture model to the DFE with a symmetric point mass of positive selection
	> This is a mixture of a 1D distribution (representing perfectly correlated fitness effects) and a 2D distribution (representing uncorrelated fitness effects)
	> The 1D distribution is a lognormal distribution (mu, sigma)
	> The 2D distribution is a bivariate lognormal distribution with shared parameters and rho = 0 (mu, sigma, rho=0)
	> The lognormal parameters mu and sigma are shared between the 1D and 2D component
	> The w parameter is the weight of the 2D distribution in the mixture model (interpreted as a DFE correlation)
	params = mu, sigma, rho, ppos, gammapos, w
	* where ppos is the proportion of positive selection
	* and gammapos is the positive selection coefficient (specified on the command line)
	"""
	print("Initializing lognormal mixture model with a symmetric point mass of positive selection...")

	# Define the DFE function used (either cache.integrate for bivariate prob. dist or cache.mixture for mixture models)
	dfe_func = DFE.mixture_symmetric_point_pos # note that the parameter order is slightly different for this case

	print("Adding a parameter for ancestral state misidentification")
	# Add a parameter for ancestral state misidentification
	dfe_func = dadi.Numerics.make_anc_state_misid_func(dfe_func)

	print("Defining selection distribution")
	# Define which selection distribution to use (a single 2D PDF for the bivariate prob. dist or a 1D and 2D PDF for the mixture model)
	dfe_dist_1d = DFE.PDFs.lognormal
	dfe_dist_2d = DFE.PDFs.biv_lognormal
	# Define the DFE dist list for the purposes of post-optimization model calculation
	dfe_dist = [dfe_dist_1d, dfe_dist_2d]

	# Define the nested and fixed indices (starting at 0)
	nested_indices = [3, 4]
	fixed_indices = [2, 4]

	# Return relevant arguments (note that dfe_dist is only returned in this way to allow automatic calculation of the model FS in the main function)
	return dfe_func, dfe_dist, nested_indices, fixed_indices

def parse_args():
	"""
	Parse the command-line arguments
	"""
	# Parse command line input
	parser = argparse.ArgumentParser()

	# Required arguments

	# This is the prefix to use
	parser.add_argument("--prefix")
	# This is the fs file name
	parser.add_argument("--fs")
	# This is the DFE model type to fit
	parser.add_argument("--dfe_model")
	# These are the names of the 1D and 2D caches
	parser.add_argument("--cache1d")
	parser.add_argument("--cache2d")
	# This is the mutation rate that we assume to convert parameter estimates
	parser.add_argument("--mut")
	# This is the effective sequence length that we assume for neutral regions (i.e., the total length of sequence from which variants *could* have been called) to convert parameter estimates
	parser.add_argument("--L_neut")
	# This is the corresponding effective sequence length that we assume for selected region
	parser.add_argument("--L_sel")
	# This is the inferred theta for the neutral case
	parser.add_argument("--theta_neut")
	
	# Required for making bootstrap replicates

	# Determine whether we need to make the bootstrapped FS
	parser.add_argument("--make_boot", default=False, action='store_true')
	# This is the VCF to use to make the bootstrap replicates
	parser.add_argument("--vcf")
	# This is the popfile to use to make the bootstrap replicates
	parser.add_argument("--popfile")

	# Required for conducting LRT

	# Determine whether we should conduct LRT
	parser.add_argument("--lrt", default=False, action='store_true')
	# This is the type of more complex model
	parser.add_argument("--complex_model")
	# These are the best-fit parameter estimates (excluding theta and the ll) for the simple DFE model
	parser.add_argument("--simple_popt", default=[], nargs="+")
	# These are the best-fit parameter estimates (excluding theta and the ll) for the complex DFE model
	parser.add_argument("--complex_popt", default=[], nargs="+")

	# Required for estimating parameter uncertainties

	# Determine whether we should estimate parameter uncertainties
	parser.add_argument("--uncert", default=False, action='store_true')
	# This specifies the complexity of the best-fit model
	parser.add_argument("--best_fit_model")
	# These are the best-fit parameter estimates (excluding theta and the ll) for the best-fit DFE model
	parser.add_argument("--best_fit_popt", default=[], nargs="+")

	# Required for CL-AIC (in addition to --best_fit_model and --best_fit_popt above)

	# Determine whether we should compute CL-AIC
	parser.add_argument("--claic", default=False, action='store_true')

	# Parse args
	args = parser.parse_args()
	return(args)

#-------------------------
# Run the main function
#-------------------------

if __name__ == '__main__':
	main()



	
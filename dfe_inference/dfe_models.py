#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# OVERVIEW
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Purpose: 
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Run as: python dfe_models.py 
# --prefix [prefix to use for output/population(s) name(s)] 
# --fs [fs file for selected sites]
# --dfe_model [argument for DFE model]
# --pneu [optional argument to fit a symmetric point mass of neutrality]
# --ppos [optional argument to fit a symmetric point mass of positive selection]
# --gammapos [specifies what gamma value the positive point mass should take]
# --cache1d [name of pickled cached spectra for 1D case] 
# --cache2d [name of pickled cached spectra for 2D case] 
# --mut [mutation rate] 
# --L_neut [effective sequence length for neutral compartment]
# --L_sel [effective sequence length for selected compartment] 
# --theta_neut [inferred theta for neutral demographic model] 
# --job_id [specifies the job number for running with CHTC]
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#--------------------------
# Import required packages
#--------------------------
import dadi
import random
import numpy as np
import matplotlib.pyplot as pyplot
import argparse
import dadi.DFE as DFE
import dill as pickle
import math

#-------------------------------------------------
# Loading the data and reading command line input
#-------------------------------------------------

def main():

	# Parse args
	args = parse_args()

	# Create variable to store the population name
	prefix = args.prefix
	# Read in the input frequency spectrum
	fs = dadi.Spectrum.from_file(args.fs)
	# Determine the DFE model type to fit
	dfe_model = args.dfe_model
	# Store a boolean that determines whether DFE models should include a point mass of neutrality
	pneu = args.pneu
	# Store a boolean that determines whether DFE models should include a point mass of positive selection
	ppos = args.ppos
	# Store the gamma value that should be used to model the point mass of positive selection
	gammapos = float(args.gammapos)
	# Load the 1D and 2D cached spectra
	cache1d = pickle.load(open(args.cache1d,'rb'))
	# print(cache1d.spectra.shape)
	# print(cache1d.spectra)
	# print(cache1d.gammas)
	cache2d = pickle.load(open(args.cache2d,'rb'))
	# print(cache2d.spectra.shape)
	# print(cache2d.spectra)
	# print(cache2d.gammas)
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
		
		print('fixed spectra:')
		print(cache1d.spectra)
		print('fixed shape:')
		print(cache1d.spectra.shape)
		print('nan entries:')
		print(np.isnan(cache1d.spectra).any(axis=(1,2)))

		print('Correcting 2D cache...')

		# Find the index of the cutoff
		cutoff_idx = min(np.where(cache2d.gammas >= correction + 1)[0].tolist())

		# Only keep indices that are less than (in abs value) the cutoff
		cache2d.spectra = cache2d.spectra[cutoff_idx:, cutoff_idx:]

		# Perform the same operation to the cache.gammas and cache.neg_gammas attributes
		cache2d.gammas = cache2d.gammas[cutoff_idx:]
		cache2d.neg_gammas = cache2d.neg_gammas[cutoff_idx:]

		print('fixed spectra:')
		print(cache2d.spectra)
		print('fixed shape:')
		print(cache2d.spectra.shape)
		print('nan entries:')
		print(np.isnan(cache2d.spectra).any(axis=(0,2,3)))
		print(np.isnan(cache2d.spectra).any(axis=(1,2,3)))

	# Create a variable to hold the specified mutation rate
	mu = float(args.mut)
	# Create a variable to hold the specified effective sequence length for neutral region
	L_neut = int(args.L_neut)
	# Create a variable to hold the specified effective sequence length for selected region
	L_sel = int(args.L_sel)
	# Parameter to hold inferred theta for the neutral case
	theta_neut = int(args.theta_neut)
	# Parameter for job ID (to be used for CHTC optimization)
	job_id = args.job_id

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

	# Construct the output file name and open it for editing
	print("Opening output file")
	if pneu:
		filename = prefix + "_" + dfe_model + "_pneu" + "_job_" + job_id + ".params"
	elif ppos:
		filename = prefix + "_" + dfe_model + "_ppos_gamma" + str(gammapos) + "_job_" + job_id + ".params"
	else:
		filename = prefix + "_" + dfe_model  + "_job_" + job_id + ".params"
	outfile = open(filename, 'w')

	print("Defining DFE function")

	# Determine which DFE model to use and collect necessary arguments for optimization function below
	# Additional nested conditionals determine whether we should fit the version of each model that includes a symmetric point mass of neutrality (controlled by the pneu argument)
	if dfe_model=="biv_lognormal_shared":
		if pneu:
			dfe_func, dfe_dist, opt_args, p0, lower_bounds, upper_bounds, fixed_params = biv_lognormal_shared_pneu(cache1d, cache2d, theta_sel)
		elif ppos:
			dfe_func, dfe_dist, opt_args, p0, lower_bounds, upper_bounds, fixed_params = biv_lognormal_shared_ppos(cache1d, cache2d, theta_sel, gammapos)
		else:
			dfe_func, dfe_dist, opt_args, p0, lower_bounds, upper_bounds, fixed_params = biv_lognormal_shared(cache1d, cache2d, theta_sel)
	elif dfe_model=="biv_lognormal_ind":
		if pneu:
			dfe_func, dfe_dist, opt_args, p0, lower_bounds, upper_bounds, fixed_params = biv_lognormal_ind_pneu(cache1d, cache2d, theta_sel)
		elif ppos:
			dfe_func, dfe_dist, opt_args, p0, lower_bounds, upper_bounds, fixed_params = biv_lognormal_ind_ppos(cache1d, cache2d, theta_sel, gammapos)
		else:
			dfe_func, dfe_dist, opt_args, p0, lower_bounds, upper_bounds, fixed_params = biv_lognormal_ind(cache1d, cache2d, theta_sel)
	elif dfe_model=="mixture_lognormal":
		if pneu:
			dfe_func, dfe_dist, opt_args, p0, lower_bounds, upper_bounds, fixed_params = mixture_lognormal_pneu(cache1d, cache2d, theta_sel)
		elif ppos:
			dfe_func, dfe_dist, opt_args, p0, lower_bounds, upper_bounds, fixed_params = mixture_lognormal_ppos(cache1d, cache2d, theta_sel, gammapos)				
		else:
			dfe_func, dfe_dist, opt_args, p0, lower_bounds, upper_bounds, fixed_params = mixture_lognormal(cache1d, cache2d, theta_sel)

	#--------------------------------
	# Optimization of DFE parameters
	#--------------------------------

	# Conduct the optimization for five replicates
	for i in range(0, 5):

		print("Conducting maximum-likelihood inference of DFE parameters (replicate " + str(i) + ")...")

		# Initialize list of likelihoods
		ll_list  = []
		# Initialize list of parameter value lists
		param_list = []

		# Obtain the best of three runs
		for j in range (0, 3):

			# Create a variable to store the likelihood
			ll_model = np.inf

			# Continue runnin the optimization with new starting values until we reach a non-Inf likelihood
			while np.isinf(ll_model):

				print("Initial starting params:")
				print(p0)

				# Preturb parameters prior to optimization
				print("Preturbing starting parameters...")
				p0 = dadi.Misc.perturb_params(p0, fold=1,
					upper_bound=upper_bounds,
					lower_bound=lower_bounds)

				print("Preturbed starting params:")
				print(p0)

				# Run optimization until we recieve a non-Inf likelihood
				print("Optimizing DFE model parameters...")
				popt = dadi.Inference.optimize(p0, fs, dfe_func, pts=None,
					func_args=opt_args,
					lower_bound=lower_bounds,
					upper_bound=upper_bounds,
					fixed_params=fixed_params,
					multinom=False)

				# Calculate the best-fit model frequency spectrum depending on the dfe_func
				print("Computing best-fit frequency spectrum...")
				if len(dfe_dist) > 1:
					# If more than one dfe_dist is specified, then provide arguments to the cache.mixture function
					model = dfe_func(popt, ns, cache1d, cache2d, dfe_dist[0], dfe_dist[1], theta_sel, pts)
				elif len(dfe_dist) == 1:
					# If only one dfe_dist is specified, then provide arguments to the cache.integrate function
					model = dfe_func(popt, ns, dfe_dist[0], theta_sel, pts)

				print("Current likelihood:")
				# Set to inf and re-run the optimization if this returns an empty likelihood
				ll_model = dadi.Inference.ll_multinom(model, fs)
				if ll_model:
					ll_model = ll_model
				else:
					ll_model = np.inf
				print(ll_model)

			# Store the likelihood
			ll_list.append(float(ll_model))
			# Store the parameter values
			param_list.append(popt)

		print(ll_list)
		print(param_list)

		# Figure out which run has the greatest fit
		max_ll = max(ll_list)
		index = ll_list.index(max_ll)

		# Use the best fit parameters as a new starting point
		p0 = param_list[index]
		print("Best params from this triplicate:")
		print(p0)
		ll_model = ll_list[index]
		print("Model likelihood...")
		print(ll_model)

		# Write results to an output file
		print("Writing results to output file...")
		# Include the p_misid in the output parameter values
		params = "\t".join(map(str, popt))
		# Append the results to the ouput file
		outfile.write(prefix + "\t" + dfe_model + "\t" + params + "\t" + str(ll_model) + "\n")

#----------------------------------
# Defining functions
#----------------------------------

def biv_lognormal_shared(cache1d, cache2d, theta_sel):
	"""
	Prepares all necessary variables to fit a bivariate lognormal distribution (with shared mu and sigma) to the DFE
	params = mu, sigma, rho
	"""
	print("Preparing for ML inference of bivariate lognormal joint DFE with shared parameters...")

	# Define the DFE function used (either cache.integrate for bivariate prob. dist or cache.mixture for mixture models)
	dfe_func = cache2d.integrate

	print("Adding a parameter for ancestral state misidentification")
	# Add a parameter for ancestral state misidentification
	dfe_func = dadi.Numerics.make_anc_state_misid_func(dfe_func)

	print("Defining selection distribution")
	# Define which selection distribution to use (a single 2D PDF for the bivariate prob. dist or a 1D and 2D PDF for the mixture model)
	dfe_dist = DFE.PDFs.biv_lognormal

	print("Defining arguments for optimization function")
	# Define arguments for the optimization function (these will depend on whether dfe_func uses cache.integrate or cache.mixture)
	opt_args = [dfe_dist, theta_sel]

	print("Adding upper, lower, and starting values for parameters (mu, sigma, rho, p_misid)")
	# Determine what the upper and lower bounds of parameter values should be and what their starting values should be
	# Ordered as [mu, sigma, rho, p_misid]
	p0 = [2, 2, 0.5, 0.05]
	print(p0)
	lower_bounds = [-15, 1e-6, 0, 0]
	print(lower_bounds)
	upper_bounds = [10, 10, 1, 0.15]
	print(upper_bounds)
	fixed_params = [None, None, None, None]
	print(fixed_params)

	# Return relevant arguments (note that dfe_dist is only returned in this way to allow automatic calculation of the model FS in the main function)
	return dfe_func, [dfe_dist], opt_args, p0, lower_bounds, upper_bounds, fixed_params

def biv_lognormal_shared_pneu(cache1d, cache2d, theta_sel):
	"""
	Prepares all necessary variables to fit a bivariate lognormal distribution (with shared mu and sigma) plus a symmetric point mass of neutrality to the DFE
	params = mu, sigma, rho, pneu, gammaneu
	* pneu defines the neutral proportion
	* gammaneu defines the neutral selection coefficient (fixed to 0)
	"""
	print("Preparing for ML inference of bivariate lognormal joint DFE with shared parameters *and* a symmetric point mass of neutrality...")

	# Define the DFE function to use (here we use the integration method for a symmetric point mass of positive selection to model a neutral component)
	dfe_func = cache2d.integrate_symmetric_point_pos

	print("Adding a parameter for ancestral state misidentification")
	# Add a parameter for ancestral state misidentification
	dfe_func = dadi.Numerics.make_anc_state_misid_func(dfe_func)

	print("Defining selection distribution")
	# Define which selection distribution to use (a single 2D PDF for the bivariate prob. dist or a 1D and 2D PDF for the mixture model)
	dfe_dist = DFE.PDFs.biv_lognormal

	print("Defining arguments for optimization function")
	# Define arguments for the optimization function (these will depend on whether dfe_func uses cache.integrate or cache.mixture)
	opt_args = [dfe_dist, theta_sel]

	print("Adding upper, lower, and starting values for parameters (mu, sigma, rho, pneu, gammaneu, p_misid)")
	# Determine what the upper and lower bounds of parameter values should be and what their starting values should be
	# Ordered as [mu, sigma, rho, pneu, gammaneu, p_misid]
	p0 = [2, 2, 0.5, 0.5, 0, 0.05]
	print(p0)
	lower_bounds = [-15, 1e-6, 0, 0, None, 0] # fixed parameters get lower bound of None
	print(lower_bounds)
	upper_bounds = [10, 10, 1, 1, None, 0.15] # fixed parameters get upper bound of None
	print(upper_bounds)
	fixed_params = [None, None, None, None, 0, None] # fixed params are gammaneu
	print(fixed_params)

	# Return relevant arguments (note that dfe_dist is only returned in this way to allow automatic calculation of the model FS in the main function)
	return dfe_func, [dfe_dist], opt_args, p0, lower_bounds, upper_bounds, fixed_params

def biv_lognormal_shared_ppos(cache1d, cache2d, theta_sel, gammapos):
	"""
	Prepares all necessary variables to fit a bivariate lognormal distribution (with shared mu and sigma) plus a symmetric point mass of positive selection to the DFE
	params = mu, sigma, rho, ppos, gammapos
	* ppos defines the positive proportion
	* gammapos defines the positive selection coefficient (fixed at the value specified on the command line)
	"""
	print("Preparing for ML inference of bivariate lognormal joint DFE with shared parameters *and* a symmetric point mass of positive selection...")

	# Define the DFE function to use (here we use the integration method for a symmetric point mass of positive selection)
	dfe_func = cache2d.integrate_symmetric_point_pos

	print("Adding a parameter for ancestral state misidentification")
	# Add a parameter for ancestral state misidentification
	dfe_func = dadi.Numerics.make_anc_state_misid_func(dfe_func)

	print("Defining selection distribution")
	# Define which selection distribution to use (a single 2D PDF for the bivariate prob. dist or a 1D and 2D PDF for the mixture model)
	dfe_dist = DFE.PDFs.biv_lognormal

	print("Defining arguments for optimization function")
	# Define arguments for the optimization function (these will depend on whether dfe_func uses cache.integrate or cache.mixture)
	opt_args = [dfe_dist, theta_sel]

	print("Adding upper, lower, and starting values for parameters (mu, sigma, rho, ppos, gammapos, p_misid)")
	# Determine what the upper and lower bounds of parameter values should be and what their starting values should be (fixing p_misid to its inferred value)
	# Ordered as [mu, sigma, rho, ppos, gammapos, p_misid]
	p0 = [2, 2, 0.5, 0.5, gammapos, 0.05]
	print(p0)
	lower_bounds = [-15, 1e-6, 0, 0, None, 0] # fixed parameters get lower bound of None
	print(lower_bounds)
	upper_bounds = [10, 10, 1, 1, None, 0.15] # fixed parameters get upper bound of None
	print(upper_bounds)
	fixed_params = [None, None, None, None, gammapos, None] # fixed params are gammapos
	print(fixed_params)

	# Return relevant arguments (note that dfe_dist is only returned in this way to allow automatic calculation of the model FS in the main function)
	return dfe_func, [dfe_dist], opt_args, p0, lower_bounds, upper_bounds, fixed_params

def biv_lognormal_ind(cache1d, cache2d, theta_sel):
	"""
	Prepares all necessary variables to fit a bivariate lognormal distribution (with distinct mu and sigma) to the DFE
	params = mu1, mu2, sigma1, sigma2, rho
	"""
	print("Preparing for ML inference of bivariate lognormal joint DFE with distinct parameters...")

	# Define the DFE function used (either cache.integrate for bivariate prob. dist or cache.mixture for mixture models)
	dfe_func = cache2d.integrate

	print("Adding a parameter for ancestral state misidentification")
	# Add a parameter for ancestral state misidentification
	dfe_func = dadi.Numerics.make_anc_state_misid_func(dfe_func)

	print("Defining selection distribution")
	# Define which selection distribution to use (a single 2D PDF for the bivariate prob. dist or a 1D and 2D PDF for the mixture model)
	dfe_dist = DFE.PDFs.biv_lognormal

	print("Defining arguments for optimization function")
	# Define arguments for the optimization function (these will depend on whether dfe_func uses cache.integrate or cache.mixture)
	opt_args = [dfe_dist, theta_sel]

	print("Adding upper, lower, and starting values for parameters (mu1, mu2, sigma1, sigma2, rho, p_misid")
	# Determine what the upper and lower bounds of parameter values should be and what their starting values should be (fixing p_misid to its inferred value)
	# Ordered as [mu1, mu2, sigma1, sigma2, rho, p_misid]
	p0 = [2, 2, 2, 2, 0.5, 0.05]
	print(p0)
	lower_bounds = [-15, -15, 1e-6, 1e-6, 0, 0]
	print(lower_bounds)
	upper_bounds = [10, 10, 10, 10, 1, 0.15]
	print(upper_bounds)
	fixed_params = [None, None, None, None, None, None]
	print(fixed_params)

	# Return relevant arguments (note that dfe_dist is only returned in this way to allow automatic calculation of the model FS in the main function)
	return dfe_func, [dfe_dist], opt_args, p0, lower_bounds, upper_bounds, fixed_params

def biv_lognormal_ind_pneu(cache1d, cache2d, theta_sel):
	"""
	Prepares all necessary variables to fit a bivariate lognormal distribution (with distinct mu and sigma) and a symmetric point mass of neutrality to the DFE
	params = mu1, mu2, sigma1, sigma2, rho, pneu, gammaneu
	* pneu defines the neutral proportion
	* gammaneu defines the neutral selection coefficient (fixed to 0)
	"""
	print("Preparing for ML inference of bivariate lognormal joint DFE with distinct parameters *and* a symmetric point mass of neutrality...")

	# Define the DFE function to use (here we use the integration method for a symmetric point mass of positive selection to model a neutral component)
	dfe_func = cache2d.integrate_symmetric_point_pos

	print("Adding a parameter for ancestral state misidentification")
	# Add a parameter for ancestral state misidentification
	dfe_func = dadi.Numerics.make_anc_state_misid_func(dfe_func)

	print("Defining selection distribution")
	# Define which selection distribution to use (a single 2D PDF for the bivariate prob. dist or a 1D and 2D PDF for the mixture model)
	dfe_dist = DFE.PDFs.biv_lognormal

	print("Defining arguments for optimization function")
	# Define arguments for the optimization function (these will depend on whether dfe_func uses cache.integrate or cache.mixture)
	opt_args = [dfe_dist, theta_sel]

	print("Adding upper, lower, and starting values for parameters (mu1, mu2, sigma1, sigma2, rho, pneu, gammaneu, p_misid")
	# Determine what the upper and lower bounds of parameter values should be and what their starting values should be (fixing p_misid to its inferred value)
	# Ordered as [mu1, mu2, sigma1, sigma2, rho, pneu, gammaneu, p_misid]
	p0 = [2, 2, 2, 2, 0.5, 0.5, 0, 0.05]
	print(p0)
	lower_bounds = [-15, -15, 1e-6, 1e-6, 0, 0, None, 0] # fixed parameters get lower bound of None
	print(lower_bounds)
	upper_bounds = [10, 10, 10, 10, 1, 1, None, 0.15] # fixed parameters get upper bound of None
	print(upper_bounds)
	fixed_params = [None, None, None, None, None, None, 0, None] # fixed parameters are gammaneu
	print(fixed_params)

	# Return relevant arguments (note that dfe_dist is only returned in this way to allow automatic calculation of the model FS in the main function)
	return dfe_func, [dfe_dist], opt_args, p0, lower_bounds, upper_bounds, fixed_params

def biv_lognormal_ind_ppos(cache1d, cache2d, theta_sel, gammapos):
	"""
	Prepares all necessary variables to fit a bivariate lognormal distribution (with distinct mu and sigma) and a symmetric point mass of positive selection to the DFE
	params = mu1, mu2, sigma1, sigma2, rho, ppos, gammapos
	* ppos defines the neutral proportion
	* gammapos defines the neutral selection coefficient (specified on the command line)
	"""
	print("Preparing for ML inference of bivariate lognormal joint DFE with distinct parameters *and* a symmetric point mass of positive selection...")

	# Define the DFE function to use (here we use the integration method for a symmetric point mass of positive selection)
	dfe_func = cache2d.integrate_symmetric_point_pos

	print("Adding a parameter for ancestral state misidentification")
	# Add a parameter for ancestral state misidentification
	dfe_func = dadi.Numerics.make_anc_state_misid_func(dfe_func)

	print("Defining selection distribution")
	# Define which selection distribution to use (a single 2D PDF for the bivariate prob. dist or a 1D and 2D PDF for the mixture model)
	dfe_dist = DFE.PDFs.biv_lognormal

	print("Defining arguments for optimization function")
	# Define arguments for the optimization function (these will depend on whether dfe_func uses cache.integrate or cache.mixture)
	opt_args = [dfe_dist, theta_sel]

	print("Adding upper, lower, and starting values for parameters (mu1, mu2, sigma1, sigma2, rho, ppos, gammapos, p_misid")
	# Determine what the upper and lower bounds of parameter values should be and what their starting values should be (fixing p_misid to its inferred value)
	# Ordered as [mu1, mu2, sigma1, sigma2, rho, ppos, gammapos, p_misid]
	p0 = [2, 2, 2, 2, 0.5, 0.5, gammapos, 0.05]
	print(p0)
	lower_bounds = [-15, -15, 1e-6, 1e-6, 0, 0, None, 0] # fixed parameters get lower bound of None
	print(lower_bounds)
	upper_bounds = [10, 10, 10, 10, 1, 1, None, 0.15] # fixed parameters get upper bound of None
	print(upper_bounds)
	fixed_params = [None, None, None, None, None, None, gammapos, None] # fixed parameters are gammapos
	print(fixed_params)

	# Return relevant arguments (note that dfe_dist is only returned in this way to allow automatic calculation of the model FS in the main function)
	return dfe_func, [dfe_dist], opt_args, p0, lower_bounds, upper_bounds, fixed_params

def mixture_lognormal(cache1d, cache2d, theta_sel):
	"""
	Prepares all necessary variables to fit a mixture model to the DFE
	> This is a mixture of a 1D distribution (representing perfectly correlated fitness effects) and a 2D distribution (representing uncorrelated fitness effects)
	> The 1D distribution is a lognormal distribution (mu, sigma)
	> The 2D distribution is a bivariate lognormal distribution with shared parameters and rho = 0 (mu, sigma, rho=0)
	> The lognormal parameters mu and sigma are shared between the 1D and 2D component
	> The w parameter is the weight of the 2D distribution in the mixture model (interpreted as a DFE correlation)
	params = mu, sigma, rho, w
	"""
	print("Preparing for ML inference of DFE parameters for a lognormal mixture model...")

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

	print("Defining arguments for optimization function")
	# Define arguments for the optimization function (these will depend on whether dfe_func uses cache.integrate or cache.mixture)
	opt_args = [cache1d, cache2d, dfe_dist_1d, dfe_dist_2d, theta_sel]

	print("Adding upper, lower, and starting values for parameters (mu, sigma, rho, w, p_misid)")
	# Determine what the upper and lower bounds of parameter values should be and what their starting values should be
	# Ordered as [mu, sigma, rho, w, p_misid]
	p0 = [2, 2, 0, 0.5, 0.05]
	print(p0)
	lower_bounds = [-15, 1e-6, None, 0, 0]
	print(lower_bounds)
	upper_bounds = [10, 10, None, 1, 0.15]
	print(upper_bounds)
	fixed_params = [None, None, 0, None, None]
	print(fixed_params)

	# Return relevant arguments (note that dfe_dist is only returned in this way to allow automatic calculation of the model FS in the main function)
	return dfe_func, dfe_dist, opt_args, p0, lower_bounds, upper_bounds, fixed_params

def mixture_lognormal_pneu(cache1d, cache2d, theta_sel):
	"""
	Prepares all necessary variables to fit a mixture model to the DFE with a symmetric point mass of neutrality
	> This is a mixture of a 1D distribution (representing perfectly correlated fitness effects) and a 2D distribution (representing uncorrelated fitness effects)
	> The 1D distribution is a lognormal distribution (mu, sigma)
	> The 2D distribution is a bivariate lognormal distribution with shared parameters and rho = 0 (mu, sigma, rho=0)
	> The lognormal parameters mu and sigma are shared between the 1D and 2D component
	> The w parameter is the weight of the 2D distribution in the mixture model (interpreted as a DFE correlation)
	params = mu, sigma, rho, pneu, gammaneu, w
	* where pneu is the proportion of neutrality
	* and gammaneu is the neutral selection coefficient (fixed at 0)
	"""
	print("Preparing for ML inference of DFE parameters for a lognormal mixture model with a symmetric point mass of neutrality...")

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

	print("Defining arguments for optimization function")
	# Define arguments for the optimization function (these will depend on whether dfe_func uses cache.integrate or cache.mixture)
	opt_args = [cache1d, cache2d, dfe_dist_1d, dfe_dist_2d, theta_sel]

	print("Adding upper, lower, and starting values for parameters (mu, sigma, rho, pneu, gammaneu, w, p_misid)")
	# Determine what the upper and lower bounds of parameter values should be and what their starting values should be (fixing p_misid to its inferred value)
	# Ordered as [mu, sigma, rho, pneu, gammaneu, w, p_misid]
	p0 = [2, 2, 0, 0.5, 0, 0.5, 0.05]
	print(p0)
	lower_bounds = [-15, 1e-6, None, 0, None, 0, 0] # fixed params have a lower bound of None
	print(lower_bounds)
	upper_bounds = [10, 10, None, 1, None, 1, 0.15] # fixed params have an upper bound of None
	print(upper_bounds)
	fixed_params = [None, None, 0, None, 0, None, None] # fixed params are 2D rho and gammaneu
	print(fixed_params)

	# Return relevant arguments (note that dfe_dist is only returned in this way to allow automatic calculation of the model FS in the main function)
	return dfe_func, dfe_dist, opt_args, p0, lower_bounds, upper_bounds, fixed_params

def mixture_lognormal_ppos(cache1d, cache2d, theta_sel, gammapos):
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
	print("Preparing for ML inference of DFE parameters for a lognormal mixture model with a symmetric point mass of positive selection...")

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

	print("Defining arguments for optimization function")
	# Define arguments for the optimization function (these will depend on whether dfe_func uses cache.integrate or cache.mixture)
	opt_args = [cache1d, cache2d, dfe_dist_1d, dfe_dist_2d, theta_sel]

	print("Adding upper, lower, and starting values for parameters (mu, sigma, rho, ppos, gammapos, w, p_misid)")
	# Determine what the upper and lower bounds of parameter values should be and what their starting values should be
	# Ordered as [mu, sigma, rho, ppos, gammapos, w, p_misid]
	p0 = [2, 2, 0, 0.5, gammapos, 0.5, 0.05]
	print(p0)
	lower_bounds = [-15, 1e-6, None, 0, None, 0, 0] # fixed params have a lower bound of None
	print(lower_bounds)
	upper_bounds = [10, 10, None, 1, None, 1, 0.15] # fixed params have an upper bound of None
	print(upper_bounds)
	fixed_params = [None, None, 0, None, gammapos, None, None] # fixed params are 2D rho and gammapos
	print(fixed_params)

	# Return relevant arguments (note that dfe_dist is only returned in this way to allow automatic calculation of the model FS in the main function)
	return dfe_func, dfe_dist, opt_args, p0, lower_bounds, upper_bounds, fixed_params

def parse_args():
	"""
	Parse the command-line arguments
	"""
	# Parse command line input
	parser = argparse.ArgumentParser()

	# This is the prefix to use
	parser.add_argument("--prefix")
	# This is the fs file name
	parser.add_argument("--fs")
	# This is the DFE model type to fit
	parser.add_argument("--dfe_model")
	# This determines whether we additionally fit a symmetric point mass of neutrality
	parser.add_argument('--pneu', default=False, action='store_true')
	# This determines whether we additionally fit a symmetric point mass of positive selection
	parser.add_argument('--ppos', default=False, action='store_true')
	# This specifies the gamma value that should be used for the positive point mass
	parser.add_argument('--gammapos', default=0)
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
	# This is the job ID for running on CHTC
	parser.add_argument("--job_id")

	# Parse args
	args = parser.parse_args()
	return(args)

#-------------------------
# Run the main function
#-------------------------

if __name__ == '__main__':
	main()



	
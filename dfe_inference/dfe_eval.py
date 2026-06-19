#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# OVERVIEW
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Purpose: the purpose of this script is to inspect the fit of DFE parameter estimates to the real data
# Output: a comparison of the empirical to model SFS and the corresponding residuals
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Run as: python dfe_eval.py 
# --prefix [prefix to use for output/population(s) name(s)] 
# --fs [fs file for selected sites]
# --dfe_model [argument for DFE model]
# --pneu [optional argument to fit a symmetric point mass of neutrality]
# --ppos [optional argument to fit a symmetric point mass of positive selection]
# --gammapos [specifies the beneficial gamma for ppos]
# --cache1d [name of pickled cached spectra for 1D case] 
# --cache2d [name of pickled cached spectra for 2D case] 
# --mut [mutation rate] 
# --L_neut [effective sequence length for neutral compartment]
# --L_sel [effective sequence length for selected compartment] 
# --theta_neut [inferred theta for neutral demographic model]
# --p_misid [inferred rate of ancestral state misidentification] 
# --popt [space-separated list of inferred parameter estimates, excluding theta and the likelihood]
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
	cache2d = pickle.load(open(args.cache2d,'rb'))
	# Placeholder to store bad gamma values
	bad_gamma_union = []

	print('Checking for NaN values in the 1D cached spectra...')
	# Check if the 1D cache contains NaN values
	if np.isnan(cache1d.spectra).any():
		print('Found NaN values (sad)...')
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
		print('Found more NaN values (double sad)...')
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
	# Parameter for ancestral state misidentification
	p_misid = float(args.p_misid)
	# Create a variable to hold the parameter estimates as a list of floats
	dfe_popt = [float(p) for p in args.popt]
	# Append the p_misid to the popt
	popt = dfe_popt + [p_misid]

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

	print("Defining DFE function")

	# Determine which DFE model to use and collect necessary arguments for optimization function below
	# Additional nested conditionals determine whether we should fit the version of each model that includes a symmetric point mass of neutrality (controlled by the pneu argument)
	if dfe_model=="biv_lognormal_shared":
		if pneu:
			dfe_func, dfe_dist, opt_args, p0, lower_bounds, upper_bounds, fixed_params = biv_lognormal_shared_pneu(cache1d, cache2d, theta_sel, p_misid)
		elif ppos:
			dfe_func, dfe_dist, opt_args, p0, lower_bounds, upper_bounds, fixed_params = biv_lognormal_shared_ppos(cache1d, cache2d, theta_sel, p_misid, gammapos)
		else:
			dfe_func, dfe_dist, opt_args, p0, lower_bounds, upper_bounds, fixed_params = biv_lognormal_shared(cache1d, cache2d, theta_sel, p_misid)
	elif dfe_model=="biv_lognormal_ind":
		if pneu:
			dfe_func, dfe_dist, opt_args, p0, lower_bounds, upper_bounds, fixed_params = biv_lognormal_ind_pneu(cache1d, cache2d, theta_sel, p_misid)
		elif ppos:
			dfe_func, dfe_dist, opt_args, p0, lower_bounds, upper_bounds, fixed_params = biv_lognormal_ind_ppos(cache1d, cache2d, theta_sel, p_misid, gammapos)
		else:
			dfe_func, dfe_dist, opt_args, p0, lower_bounds, upper_bounds, fixed_params = biv_lognormal_ind(cache1d, cache2d, theta_sel, p_misid)
	elif dfe_model=="mixture_lognormal":
		if pneu:
			dfe_func, dfe_dist, opt_args, p0, lower_bounds, upper_bounds, fixed_params = mixture_lognormal_pneu(cache1d, cache2d, theta_sel, p_misid)
		elif ppos:
			dfe_func, dfe_dist, opt_args, p0, lower_bounds, upper_bounds, fixed_params = mixture_lognormal_ppos(cache1d, cache2d, theta_sel, p_misid, gammapos)				
		else:
			dfe_func, dfe_dist, opt_args, p0, lower_bounds, upper_bounds, fixed_params = mixture_lognormal(cache1d, cache2d, theta_sel, p_misid)

	#--------------------------------
	# Model evaluation
	#--------------------------------

	print("Best-fit DFE parameter values...")
	print(popt)

	# For the cases where we're fitting one of the non-neutral DFE models
	if dfe_model!='neutral':
		# Calculate the best-fit model frequency spectrum depending on the dfe_func
		print("Computing model frequency spectrum for specified DFE parameters...")
		if len(dfe_dist) > 1:
			# If more than one dfe_dist is specified, then provide arguments to the cache.mixture function
			model = dfe_func(popt, ns, cache1d, cache2d, dfe_dist[0], dfe_dist[1], theta_sel, pts)
		elif len(dfe_dist) == 1:
			# If only one dfe_dist is specified, then provide arguments to the cache.integrate function
			model = dfe_func(popt, ns, dfe_dist[0], theta_sel, pts)
		#print(model)

	# For the case of 100% neutrality
	elif dfe_model=='neutral':
		# Calculate the best-fit model frequency spectrum under neutrality
		print("Computing model frequency spectrum for fully neutral DFE...")
		# To do this, we get the slice of the cache that corresponds to a gamma=0 (see the API for the pos selection integrate functions for reference)
		exp_sfs = cache2d.spectra[cache2d.gammas == 0, cache2d.gammas == 0][0]
		# To get the model expectations, multiply this FS by the selected theta
		model = exp_sfs*theta_sel
		# Then model ancestral state misidentification
		model = dadi.Numerics.apply_anc_state_misid(model, p_misid)
		# And turn it into a dadi spectrum object
		model = dadi.Spectrum(model)
		#print(model)		

	print("Model likelihood...")
	ll_model = dadi.Inference.ll_multinom(model, fs)
	print(ll_model)

	#------------------------------------
	# Generating residuals and DFE plots
	#------------------------------------

	print("Making residual plots...")

	# Open the figure
	pylab.figure(figsize=(7,6))

	# Assign vmin and vmax for plotting 2D
	vmax=np.max(np.maximum(model, fs))
	print(vmax)
	#vmin=np.min(np.minimum(model, fs))
	vmin=1 # Apparently this will throw and error if zero
	print(vmin)

	# Plot the 2D comparison
	dadi.Plotting.plot_2d_comp_multinom(model, fs, vmin=vmin, vmax=vmax, show=False, resid_range=5)
	# Save the figure
	if pneu:
		pylab.savefig((str(prefix) + "_" + str(dfe_model) + "_pneu" + "_residuals.png"), dpi=350)
	elif ppos:
		pylab.savefig((str(prefix) + "_" + str(dfe_model) + "_ppos_gamma" + str(gammapos) + "_residuals.png"), dpi=350)
	else:
		pylab.savefig((str(prefix) + "_" + str(dfe_model) + "_residuals.png"), dpi=350)

	# Open another figure
	pylab.figure(figsize=(7,6))
	# Marginalize over the second population
	data_fs1 = fs.marginalize([1])
	model_fs1 = model.marginalize([1])
	# Plot the 1D comparison between model and data
	plot_1d_comp_multinom(model_fs1, data_fs1, "Population 1", dfe_model, data_fs1.sample_sizes, show=False)
	# Save the figure
	if pneu:
		pylab.savefig((str(prefix) + "_" + str(dfe_model) + "_pneu" + "_pop1_1Dresiduals.png"), dpi=350)
	elif ppos:
		pylab.savefig((str(prefix) + "_" + str(dfe_model) + "_ppos_gamma" + str(gammapos) + "_pop1_1Dresiduals.png"), dpi=350)
	else:
		pylab.savefig((str(prefix) + "_" + str(dfe_model) + "_pop1_1Dresiduals.png"), dpi=350)

	# Open another figure
	pylab.figure(figsize=(7,6))
	# Marginalize over the first population
	data_fs2 = fs.marginalize([0])
	model_fs2 = model.marginalize([0])
	# Plot the 1D comparison between model and data
	plot_1d_comp_multinom(model_fs2, data_fs2, "Population 2", dfe_model, data_fs2.sample_sizes, show=False)
	# Save the figure
	if pneu:
		pylab.savefig((str(prefix) + "_" + str(dfe_model) + "_pneu" + "_pop2_1Dresiduals.png"), dpi=350)
	elif ppos:
		pylab.savefig((str(prefix) + "_" + str(dfe_model) + "_ppos_gamma" + str(gammapos) + "_pop2_1Dresiduals.png"), dpi=350)
	else:
		pylab.savefig((str(prefix) + "_" + str(dfe_model) + "_pop2_1Dresiduals.png"), dpi=350)

#----------------------------------
# Defining functions
#----------------------------------
def biv_lognormal_shared(cache1d, cache2d, theta_sel, p_misid):
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

	print("Adding upper, lower, and starting values for parameters (mu, sigma, rho, p_misid) (IGNORED FOR MODEL EVAL)")
	# Determine what the upper and lower bounds of parameter values should be and what their starting values should be (fixing p_misid to its inferred value)
	# Ordered as [mu, sigma, rho, pneu, gammaneu, p_misid]
	p0 = [2, 2, 0.5, p_misid]
	print(p0)
	lower_bounds = [1e-6, 1e-6, 0, None]
	print(lower_bounds)
	upper_bounds = [10, 5, 1, None]
	print(upper_bounds)
	fixed_params = [None, None, None, p_misid]
	print(fixed_params)

	# Return relevant arguments (note that dfe_dist is only returned in this way to allow automatic calculation of the model FS in the main function)
	return dfe_func, [dfe_dist], opt_args, p0, lower_bounds, upper_bounds, fixed_params

def biv_lognormal_shared_pneu(cache1d, cache2d, theta_sel, p_misid):
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

	print("Adding upper, lower, and starting values for parameters (mu, sigma, rho, pneu, gammaneu, p_misid) (IGNORED FOR MODEL EVAL)")
	# Determine what the upper and lower bounds of parameter values should be and what their starting values should be (fixing p_misid to its inferred value)
	# Ordered as [mu, sigma, rho, pneu, gammaneu, p_misid]
	p0 = [2, 2, 0.5, 0.5, 0, p_misid]
	print(p0)
	lower_bounds = [1e-6, 1e-6, 0, 0, None, None] # fixed parameters get lower bound of None
	print(lower_bounds)
	upper_bounds = [10, 5, 1, 1, None, None] # fixed parameters get upper bound of None
	print(upper_bounds)
	fixed_params = [None, None, None, None, 0, p_misid] # fixed params are gammaneu and p_misid
	print(fixed_params)

	# Return relevant arguments (note that dfe_dist is only returned in this way to allow automatic calculation of the model FS in the main function)
	return dfe_func, [dfe_dist], opt_args, p0, lower_bounds, upper_bounds, fixed_params

def biv_lognormal_shared_ppos(cache1d, cache2d, theta_sel, p_misid, gammapos):
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
	p0 = [2, 2, 0.5, 0.5, gammapos, p_misid]
	print(p0)
	lower_bounds = [1e-6, 1e-6, 0, 0, None, None] # fixed parameters get lower bound of None
	print(lower_bounds)
	upper_bounds = [10, 5, 1, 1, None, None] # fixed parameters get upper bound of None
	print(upper_bounds)
	fixed_params = [None, None, None, None, gammapos, p_misid] # fixed params are gammapos and p_misid
	print(fixed_params)

	# Return relevant arguments (note that dfe_dist is only returned in this way to allow automatic calculation of the model FS in the main function)
	return dfe_func, [dfe_dist], opt_args, p0, lower_bounds, upper_bounds, fixed_params

def biv_lognormal_ind(cache1d, cache2d, theta_sel, p_misid):
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

	print("Adding upper, lower, and starting values for parameters (mu1, mu2, sigma1, sigma2, rho, p_misid) (IGNORED FOR MODEL EVAL)")
	# Determine what the upper and lower bounds of parameter values should be and what their starting values should be (fixing p_misid to its inferred value)
	# Ordered as [mu1, mu2, sigma1, sigma2, rho, p_misid]
	p0 = [2, 2, 2, 2, 0.5, p_misid]
	print(p0)
	lower_bounds = [1e-6, 1e-6, 1e-6, 1e-6, 0, None]
	print(lower_bounds)
	upper_bounds = [10, 10, 5, 5, 1, None]
	print(upper_bounds)
	fixed_params = [None, None, None, None, None, p_misid]
	print(fixed_params)

	# Return relevant arguments (note that dfe_dist is only returned in this way to allow automatic calculation of the model FS in the main function)
	return dfe_func, [dfe_dist], opt_args, p0, lower_bounds, upper_bounds, fixed_params

def biv_lognormal_ind_pneu(cache1d, cache2d, theta_sel, p_misid):
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

	print("Adding upper, lower, and starting values for parameters (mu1, mu2, sigma1, sigma2, rho, pneu, gammaneu, p_misid) (IGNORED FOR MODEL EVAL)")
	# Determine what the upper and lower bounds of parameter values should be and what their starting values should be (fixing p_misid to its inferred value)
	# Ordered as [mu1, mu2, sigma1, sigma2, rho, pneu, gammaneu, p_misid]
	p0 = [2, 2, 2, 2, 0.5, 0.5, 0, p_misid]
	print(p0)
	lower_bounds = [1e-6, 1e-6, 1e-6, 1e-6, 0, 0, None, None] # fixed parameters get lower bound of None
	print(lower_bounds)
	upper_bounds = [10, 10, 5, 5, 1, 1, None, None] # fixed parameters get upper bound of None
	print(upper_bounds)
	fixed_params = [None, None, None, None, None, None, 0, p_misid] # fixed parameters are gammaneu and p_misid
	print(fixed_params)

	# Return relevant arguments (note that dfe_dist is only returned in this way to allow automatic calculation of the model FS in the main function)
	return dfe_func, [dfe_dist], opt_args, p0, lower_bounds, upper_bounds, fixed_params

def biv_lognormal_ind_ppos(cache1d, cache2d, theta_sel, p_misid, gammapos):
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
	p0 = [2, 2, 2, 2, 0.5, 0.5, gammapos, p_misid]
	print(p0)
	lower_bounds = [1e-6, 1e-6, 1e-6, 1e-6, 0, 0, None, None] # fixed parameters get lower bound of None
	print(lower_bounds)
	upper_bounds = [10, 10, 5, 5, 1, 1, None, None] # fixed parameters get upper bound of None
	print(upper_bounds)
	fixed_params = [None, None, None, None, None, None, gammapos, p_misid] # fixed parameters are gammapos and p_misid
	print(fixed_params)

	# Return relevant arguments (note that dfe_dist is only returned in this way to allow automatic calculation of the model FS in the main function)
	return dfe_func, [dfe_dist], opt_args, p0, lower_bounds, upper_bounds, fixed_params

def mixture_lognormal(cache1d, cache2d, theta_sel, p_misid):
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

	print("Adding upper, lower, and starting values for parameters (mu, sigma, rho, w, p_misid) (IGNORED FOR MODEL EVAL)")
	# Determine what the upper and lower bounds of parameter values should be and what their starting values should be (fixing p_misid to its inferred value)
	# Ordered as [mu, sigma, rho, w, p_misid]
	p0 = [2, 2, 0, 0.5, p_misid]
	print(p0)
	lower_bounds = [1e-6, 1e-6, None, 0, None]
	print(lower_bounds)
	upper_bounds = [10, 5, None, 1, None]
	print(upper_bounds)
	fixed_params = [None, None, 0, None, p_misid]
	print(fixed_params)

	# Return relevant arguments (note that dfe_dist is only returned in this way to allow automatic calculation of the model FS in the main function)
	return dfe_func, dfe_dist, opt_args, p0, lower_bounds, upper_bounds, fixed_params

def mixture_lognormal_pneu(cache1d, cache2d, theta_sel, p_misid):
	"""
	Prepares all necessary variables to fit a mixture model to the DFE with a symmetric point mass of neutrality
	> This is a mixture of a 1D distribution (representing perfectly correlated fitness effects) and a 2D distribution (representing uncorrelated fitness effects)
	> The 1D distribution is a lognormal distribution (mu, sigma)
	> The 2D distribution is a bivariate lognormal distribution with shared parameters and rho = 0 (mu, sigma, rho=0)
	> The lognormal parameters mu and sigma are shared between the 1D and 2D component
	> The w parameter is the weight of the 2D distribution in the mixture model (interpreted as a DFE correlation)
	params = mu, sigma, rho, pneu, gammaneu, w
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

	print("Adding upper, lower, and starting values for parameters (mu, sigma, rho, pneu, gammaneu, w, p_misid) (IGNORED FOR MODEL EVAL)")
	# Determine what the upper and lower bounds of parameter values should be and what their starting values should be (fixing p_misid to its inferred value)
	# Ordered as [mu, sigma, rho, pneu, gammaneu, w, p_misid]
	p0 = [2, 2, 0, 0.5, 0, 0.5, p_misid]
	print(p0)
	lower_bounds = [1e-6, 1e-6, None, 0, None, 0, None] # fixed params have a lower bound of None
	print(lower_bounds)
	upper_bounds = [10, 5, None, 1, None, 1, None] # fixed params have an upper bound of None
	print(upper_bounds)
	fixed_params = [None, None, 0, None, 0, None, p_misid] # fixed params are 2D rho, gammaneu, and p_misid
	print(fixed_params)

	# Return relevant arguments (note that dfe_dist is only returned in this way to allow automatic calculation of the model FS in the main function)
	return dfe_func, dfe_dist, opt_args, p0, lower_bounds, upper_bounds, fixed_params

def mixture_lognormal_ppos(cache1d, cache2d, theta_sel, p_misid, gammapos):
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
	# Determine what the upper and lower bounds of parameter values should be and what their starting values should be (fixing p_misid to its inferred value)
	# Ordered as [mu, sigma, rho, ppos, gammapos, w, p_misid]
	p0 = [2, 2, 0, 0.5, gammapos, 0.5, p_misid]
	print(p0)
	lower_bounds = [1e-6, 1e-6, None, 0, None, 0, None] # fixed params have a lower bound of None
	print(lower_bounds)
	upper_bounds = [10, 5, None, 1, None, 1, None] # fixed params have an upper bound of None
	print(upper_bounds)
	fixed_params = [None, None, 0, None, gammapos, None, p_misid] # fixed params are 2D rho, gammapos, and p_misid
	print(fixed_params)

	# Return relevant arguments (note that dfe_dist is only returned in this way to allow automatic calculation of the model FS in the main function)
	return dfe_func, dfe_dist, opt_args, p0, lower_bounds, upper_bounds, fixed_params

def plot_1d_comp_Poisson(model, data, pop_name, model_id, ns, folded, fig_num=None, residual='Anscombe', plot_masked=False, show=True):
	"""
	Poisson comparison between 1d model and data.


	model: 1-dimensional model SFS
	data: 1-dimensional data SFS
	fig_num: Clear and use figure fig_num for display. If None, an new figure
			 window is created.
	residual: 'Anscombe' for Anscombe residuals, which are more normally
			  distributed for Poisson sampling. 'linear' for the linear
			  residuals, which can be less biased.
	plot_masked: Additionally plots (in open circles) results for points in the 
				 model or data that were masked.
	show: If True, execute pylab.show command to make sure plot displays.
	"""
	if fig_num is None:
		f = pylab.gcf()
	else:
		f = pylab.figure(fig_num, figsize=(7,7))
	pylab.clf()

	if data.folded and not model.folded:
		model = model.fold()

	masked_model, masked_data = dadi.Numerics.intersect_masks(model, data)

	ax = pylab.subplot(2,1,1)
	pylab.semilogy(masked_data, '-ob', label='data')
	pylab.semilogy(masked_model, '-or', label='model')

	if plot_masked:
		pylab.semilogy(masked_data.data, '--ob', mfc='w', zorder=-100)
		pylab.semilogy(masked_model.data, '--or', mfc='w', zorder=-100)

	ax.legend(loc='lower left')

	res_ax = pylab.subplot(2,1,2, sharex = ax)
	if residual == 'Anscombe':
		resid = dadi.Inference.Anscombe_Poisson_residual(masked_model, masked_data)
	elif residual == 'linear':
		resid = dadi.Inference.linear_Poisson_residual(masked_model, masked_data)
	else:
		raise ValueError("Unknown class of residual '%s'." % residual)
	pylab.plot(resid, '-og')
	if plot_masked:
		pylab.plot(resid.data, '--og', mfc='w', zorder=-100)

	# Edited portion (10/10/25)-------------------------------------------------------------------
	if (folded):
		print('Using folded 1D spectrum for plotting...')
		ax.set_xlim(0, data.shape[0]-(ns/2))
	else: 
		print('Using unfolded 1D spectrum for plotting...')
		ax.set_xlim(0, data.shape[0]) 

	f.suptitle((str(pop_name) + "; " + str(model_id) + " model fit"))
	ax.set_xlabel('Derived allele count')
	ax.set_ylabel('Number of SNPs')
	res_ax.set_xlabel('Derived allele count')
	res_ax.set_ylabel('Residuals')

	pyplot.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.5)
	#----------------------------------------------------------------------------------------------

	if show:
		pylab.show()

def plot_1d_comp_multinom(model, data, pop_name, model_id, ns, folded=False, fig_num=None, residual='Anscombe', plot_masked=False, show=True):
	"""
	Mulitnomial comparison between 1d model and data.


	model: 1-dimensional model SFS
	data: 1-dimensional data SFS
	fig_num: Clear and use figure fig_num for display. If None, an new figure
			 window is created.
	residual: 'Anscombe' for Anscombe residuals, which are more normally
			  distributed for Poisson sampling. 'linear' for the linear
			  residuals, which can be less biased.
	plot_masked: Additionally plots (in open circles) results for points in the 
				 model or data that were masked.
	show: If True, execute pylab.show command to make sure plot displays.

	This comparison is multinomial in that it rescales the model to optimally
	fit the data.
	"""
	model = dadi.Inference.optimally_scaled_sfs(model, data)

	plot_1d_comp_Poisson(model, data, pop_name, model_id, ns, folded, fig_num, residual,
						 plot_masked, show)

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
	# This is the inferred rate of ancestral state misidentification
	parser.add_argument("--p_misid")
	# This is the job ID for running on CHTC
	parser.add_argument("--job_id")
	# These are the best-fit parameter estimates (excluding theta and the ll) for the given model
	parser.add_argument("--popt", default=[], nargs="+")

	# Parse args
	args = parser.parse_args()
	return(args)

#-------------------------
# Run the main function
#-------------------------

if __name__ == '__main__':
	main()



	
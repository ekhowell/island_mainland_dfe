#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# OVERVIEW
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Purpose: the purpose of this script is to inspect the fit of demographic parameter estimates to the real data
# Output: a comparison of the empirical to model SFS and the corresponding residuals
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Run as: python demog_eval.py --prefix [population name(s)] --fs [fs file] --model [model name] --popt [param1 param2 ... paramN] --p_misid
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#--------------------------
# Import required packages
#--------------------------
import dadi
import random
import numpy as np
import matplotlib.pyplot as pyplot
import argparse
import re
import pylab

#----------------------------------
# DEFINE THE MAIN SCRIPT FUNCTION
#----------------------------------

def main():

	# Parse command-line arguments
	args = parse_args()

	# Create variable to store the population name
	prefix = args.prefix
	# Read in the input frequency spectrum
	fs = dadi.Spectrum.from_file(args.fs)
	# Create a variable to hold the model ID
	model_id = args.model
	# Create a variable to hold the parameter estimates
	popt_est = args.popt
	# Create a variable to hold the parameter estimates as a list of floats
	popt = [float(p) for p in popt_est]
	# Parameter for ancestral state misidentification
	p_misid = args.p_misid

	# Print out the best-fit model parameters
	print("Best-fit " + str(prefix) + " parameters for the " + str(model_id) + " model: " + ", ".join(map(str, popt_est)))

	# Determine what the sample size is based on the input fs
	ns = fs.sample_sizes

	# Need to define the grid points for optimization
	pts = [x for x in range(0, len(ns)*10, 10)] + (max(ns)*2)

	# Need to determine pop_ids based on input fs
	pop_ids = fs.pop_ids

	# Wrap our demographic function based on user specified model
	if (model_id=="boston_harbor"):
		model_func = boston_harbor
	elif (model_id=="gulf_islands_ivi"):
		model_func = gulf_islands_ivi
	elif (model_id=="gulf_islands_ivm"):
		model_func = gulf_islands_ivm

	# Determine whether to add a parameter for ancestral state misidentification
	if p_misid:
		print("Including a parameter for ancestral state misidentification")
		# If so, add this parameter to the model
		model_func = dadi.Numerics.make_anc_state_misid_func(model_func)

	# Wrap the demographic model in a function that utilizes grid points which increases dadi's ability to more accurately generate a model frequency spectrum.
	demo_model_ex = dadi.Numerics.make_extrap_func(model_func)

	print("Computing best-fit frequency spectrum...")

	# Calculate the best-fit model FS
	model = demo_model_ex(popt, ns, pts)

	# Calculate the likelihood of the data given the model FS
	ll_model = dadi.Inference.ll_multinom(model, fs)
				
	# Calculate the corresponding theta value
	theta = dadi.Inference.optimal_sfs_scaling(model, fs)

	# Print out the computed likelihood
	print("Maximum composite likelihood: " + str(ll_model))

	# Print out the optimal value of theta
	print("Optimal value of theta: " + str(theta))

	# Plot a comparison of the model FS with the data

	# Open the figure
	pylab.figure(figsize=(7,6))

	# Plot the 2D comparison
	dadi.Plotting.plot_2d_comp_multinom(model, fs, vmin=1, vmax=fs.max(), show=False)
	# Save the figure
	pylab.savefig((str(prefix) + "_" + str(model_id) + "_residuals.png"), dpi=250)

	# Open another figure
	pylab.figure(figsize=(7,6))
	# Marginalize over the second population
	data_fs1 = fs.marginalize([1])
	model_fs1 = model.marginalize([1])
	# Plot the 1D comparison between model and data
	plot_1d_comp_multinom(model_fs1, data_fs1, "Population 1", model_id, data_fs1.sample_sizes, folded=False, show=False)
	# Save the figure
	pylab.savefig((str(prefix) + "_" + str(model_id) + "_pop1_1Dresiduals.png"), dpi=250)

	# Open another figure
	pylab.figure(figsize=(7,6))
	# Marginalize over the first population
	data_fs2 = fs.marginalize([0])
	model_fs2 = model.marginalize([0])
	# Plot the 1D comparison between model and data
	plot_1d_comp_multinom(model_fs2, data_fs2, "Population 2", model_id, data_fs2.sample_sizes, folded=False, show=False)
	# Save the figure
	pylab.savefig((str(prefix) + "_" + str(model_id) + "_pop2_1Dresiduals.png"), dpi=250)

#---------------------
# DEFINE FUNCTIONS
#---------------------

# Best-fit demographic model for Boston Harbor in dadi-format
def boston_harbor(params, ns, pts):
	"""
	Parameter values:
	1. nu1 = the ratio of population 1 size to the ancestral population size
	2. nu2 = the ratio of population 2 size to the ancestral population size
	3. T = the time since the population split
	Overview: the ancestral population splits into population 1 with size nu1 and population 2 with size nu2 T generations from the present.
	"""
	# Initialize the parameter values
	nu1, nu2, T = params	
	# Specify the density of the grid over which phi is evaluated
	xx = dadi.Numerics.default_grid(pts)
	# Ancestral dynamics #
	# Create a initial 1D phi that represents our equilibrium ancestral population with a relative population size of 1
	phi = dadi.PhiManip.phi_1D(xx)
	# Epoch 1 dynamics #
	# Split the ancestral population into two nascent populations
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	# After the splitting, the populations take constant sizes nu1 and nu2 for T generations
	phi = dadi.Integration.two_pops(phi, xx, T, nu1, nu2)
	# Return the frequency spectrum obtained from the model
	fs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))
	return fs

# Define the demographic function for a population split followed by constant sizes for both population 1 and 2 with symmetric migration
# This corresponds to the Saturna–Pender comparison
def gulf_islands_ivi(params, ns, pts):
	"""
	Parameter values:
	1. nu1 = the ratio of population 1's size to the ancestral population size
	2. nu2 = the ratio of population 2's size to the ancestral population size
	3. T = the time since the population split
	4. m = the symmetric migration rate
	Overview: the ancestral population splits into population 1 with size nu1 and population 2 with size nu2 T generations from the present.
		After the split, there is symmetric migration at rate m.
	[NOTE] all parameters are given relative to the reference size. The ancestral reference size must be calculated from theta
	"""
	# Initialize the parameter values
	nu1, nu2, T, m = params
	# Specify the density of the grid over which phi is evaluated
	xx = dadi.Numerics.default_grid(pts)
	# Ancestral dynamics #
	# Create a initial 1D phi that represents our equilibrium ancestral population with a relative population size of 1
	phi = dadi.PhiManip.phi_1D(xx)
	# Epoch 1 dynamics #
	# Split the ancestral population into two nascent populations 
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	# After the splitting, the populations take constant sizes nu1 and nu2 for T generations with symmetric migration rate m
	phi = dadi.Integration.two_pops(phi, xx, T, nu1, nu2, m12=m, m21=m)
	# Return the frequency spectrum obtained from the model
	fs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))
	return fs

# Define the demographic function for a population split followed by a post-split size change with multiple mig parameters
# This corresponds to the Saturna–Maple Ridge and Pender–Maple Ridge comparisons
def gulf_islands_ivm(params, ns, pts):
	"""
	Parameter values:
	1. nu1E1 = the ratio of population 1 size to the ancestral population size during first epoch
	2. nu2E1 = the ratio of population 2 size to the ancestral population size during the first epoch
	3. nu1E2 = the ratio of population 1 size to the ancestral population size during second epoch
	4. nu2E2 = the ratio of population 2 size to the ancestral population size during the second epoch
	5. TE1 = the duration of the first epoch
	6. TE2 = the duration of the second epoch
	7. mE1 = the migration rate during the first epoch
	8. mE2 = the migration rate during the second epoch
	Overview: the ancestral population splits into population 1 with size nu1E1 and population 2 with size nu2E1 and migration rate mE1. After TE1 generations, the populations enter a second epoch with sizes nu1E1 and nu2E2 that lasts for TE2 generations with migration rate mE2.
	[NOTE] all parameters are given relative to the reference size. The ancestral reference size must be calculated from theta
	"""
	# Initialize the parameter values
	nu1E1, nu2E1, nu1E2, nu2E2, TE1, TE2, mE1, mE2 = params
	# Specify the density of the grid over which phi is evaluated
	xx = dadi.Numerics.default_grid(pts)
	# Ancestral dynamics #
	# Create a initial 1D phi that represents our equilibrium ancestral population with a relative population size of 1
	phi = dadi.PhiManip.phi_1D(xx)
	# Epoch 1 dynamics #
	# Split the ancestral population into two nascent populations that represent the mainland ancestor and the island ancestor
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	# After the splitting, the populations take constant sizes nu1E1 and nu2E1 for TE1 generations with symmetric migration rate mE1
	phi = dadi.Integration.two_pops(phi, xx, TE1, nu1E1, nu2E1, m12=mE1, m21=mE1)
	# Epoch 2 dynamics #
	# The populations then take constant sizes nu1E2 and nu2E2 for TE2 generations with symmetric migration rate mE2
	phi = dadi.Integration.two_pops(phi, xx, TE2, nu1E2, nu2E2, m12=mE2, m21=mE2)
	# Return the frequency spectrum obtained from the model
	fs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))
	return fs

def plot_1d_comp_Poisson(model, data, prefix, model_id, ns, folded, fig_num=None, residual='Anscombe', plot_masked=False, show=True):
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

	f.suptitle((str(prefix) + "; " + str(model_id) + " model fit"))
	ax.set_xlabel('Derived allele count')
	ax.set_ylabel('Number of SNPs')
	res_ax.set_xlabel('Derived allele count')
	res_ax.set_ylabel('Residuals')

	pyplot.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.5)
	#----------------------------------------------------------------------------------------------

	if show:
		pylab.show()

def plot_1d_comp_multinom(model, data, prefix, model_id, ns, folded, fig_num=None, residual='Anscombe', plot_masked=False, show=True):
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

	plot_1d_comp_Poisson(model, data, prefix, model_id, ns, folded, fig_num, residual,
						 plot_masked, show)

def parse_args():
	"""
	Parse the command-line arguments
	"""
	# Parse command line input
	parser = argparse.ArgumentParser()

	# Prefix to use
	parser.add_argument("--prefix")
	# This is the fs file name
	parser.add_argument("--fs")
	# This is the demographic model we wish to fit 
	parser.add_argument("--model")
	# These are the best-fit parameter estimates (excluding theta and the likelihood) for the given model
	parser.add_argument("--popt", default=[], nargs="+")
	# Indicates whether to add a parameter for ancestral state misidentification
	parser.add_argument("--p_misid", action='store_true', default=False)

	# Parse args
	args = parser.parse_args()
	return(args)

#-------------------------
# RUN THE MAIN FUNCTION
#-------------------------

if __name__ == '__main__':
	main()








	
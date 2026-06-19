#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# OVERVIEW
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Purpose: the purpose of this script is to obtain the optimal parameter estimates for the given demograhic model
# Output: a tab-delimited file with the rows formatted as follows: [Model] [Param 1] [Param 2] ... [Theta] [Likelihood] [Popt]
# Where [Popt] represents the un-converted parameter estimates
# Convergence and best-fit parameter values can then be evaluated in R
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Run as: python demog_models.py --prefix [prefix to use for output/population(s) name(s)] --fs [fs file] --model [model name] --opt_num [optimization trial number] --mut [mutation rate] --L [effective sequence length] --p_misid [include a parameter for ancestral state misidentification?] 
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#--------------------------
# Import required packages
#--------------------------
import dadi
import random
import numpy as np
import matplotlib.pyplot as pyplot
import argparse

#---------------------------------------------
# Loading the data and initializing variables
#---------------------------------------------

def main():

	# Parse args
	args = parse_args()

	# Create variable to store the population name
	prefix = args.prefix
	# Read in the input frequency spectrum
	fs = dadi.Spectrum.from_file(args.fs)
	# Create a variable to hold the model ID
	model_id = args.model
	# Create a variable to hold which independent optimization trial this is
	opt_num = args.opt_num
	# Create a variable to hold the specified mutation rate
	mu = float(args.mut)
	# Create a variable to hold the specified effective sequence length
	seqL = int(args.L)
	# Parameter for ancestral state misidentification
	p_misid = args.p_misid

	#------------------------
	# Simulation and fitting
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

	# Construct the output file name and open it for editing
	if p_misid:
		filename = prefix + "_" + model_id  + "_ancMisid" + "_optRun" + opt_num + ".params"
	else:
		filename = prefix + "_" + model_id + "_optRun" + opt_num + ".params"
	outfile = open(filename, 'w')

	# Wrap our demographic function based on user specified model
	if model_id=="boston_harbor":
		model_func = boston_harbor
	elif (model_id=="gulf_islands_ivi"):
		model_func = gulf_islands_ivi
	elif (model_id=="gulf_islands_ivm"):
		model_func = gulf_islands_ivm

	# Determine whether to add a parameter for ancestral state misidentification
	if p_misid:
		print("Adding a parameter for ancestral state misidentification")
		# If so, add this parameter to the model
		model_func = dadi.Numerics.make_anc_state_misid_func(model_func)

	# Wrap the demographic model in a function that utilizes grid points which increases dadi's ability to more accurately generate a model frequency spectrum.
	demo_model_ex = dadi.Numerics.make_extrap_func(model_func)

	# Determine what the upper and lower bounds of parameter values should be and what their starting values should be
	if model_id=="boston_harbor":
		upper_bound,lower_bound,p0 = instantiate_boston_harbor()
	elif (model_id=="gulf_islands_ivi"):
		upper_bound,lower_bound,p0 = instantiate_gulf_islands_ivi()
	elif (model_id=="gulf_islands_ivm"):
		upper_bound,lower_bound,p0 = instantiate_gulf_islands_ivm()

	# Determine whether we are modeling ancestral state misidentification
	if p_misid:
		print("Adding upper, lower, and starting values for p_misid")
		# If so, add this parameter to the p0 starting parameters
		upper_bound.append(1)
		lower_bound.append(0)
		p0.append(1e-4)

	#------------------------
	# Optimization routine
	#------------------------
	print("Conducting maximum-likelihood inference of parameters...")

	# Preturb parameters prior to optimization
	print("Preturbing starting parameters...")
	p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

	# Conduct MLE 
	print("Optimizing demographic parameters...")
	popt, ll_model = dadi.Inference.opt(p0, fs, demo_model_ex, pts, lower_bound=lower_bound, upper_bound=upper_bound)

	# Calculate the best-fit model frequency spectrum
	print("Computing best-fit frequency spectrum...")
	model = demo_model_ex(popt, ns, pts)
	# Additionally, return the optimal value of theta given the model
	theta = dadi.Inference.optimal_sfs_scaling(model, fs)

	# I'm going to do the parameter conversion here for simplicity (requires mu and effective sequence length)
	# Determine how the parameter conversion should be performed (depends on the model)
	if p_misid:
		if (model_id=="boston_harbor"):
			converted = params_boston_harbor(theta, mu, seqL, popt) + [float(popt[-1])]
			param_vals = "\t".join(map(str, converted))
		elif (model_id=="gulf_islands_ivi"):
			converted = params_gulf_islands_ivi(theta, mu, seqL, popt) + [float(popt[-1])]
			param_vals = "\t".join(map(str, converted))
		elif (model_id=="gulf_islands_ivm"):
			converted = params_gulf_islands_ivm(theta, mu, seqL, popt) + [float(popt[-1])]
			param_vals = "\t".join(map(str, converted))
	else:
		if (model_id=="boston_harbor"):
			converted = params_boston_harbor(theta, mu, seqL, popt)
			param_vals = "\t".join(map(str, converted))
		elif (model_id=="gulf_islands_ivi"):
			converted = params_gulf_islands_ivi(theta, mu, seqL, popt)
			param_vals = "\t".join(map(str, converted))
		elif (model_id=="gulf_islands_ivm"):
			converted = params_gulf_islands_ivm(theta, mu, seqL, popt)
			param_vals = "\t".join(map(str, converted))

	# Append the results to the ouput file
	outfile.write(prefix + "\t" + model_id + "\t" + param_vals + "\t" + str(ll_model) + "\t" + "\t".join(map(str, popt)) + "\n")
	outfile.close()


#----------------------------------
# Specifying the demographic model
#----------------------------------

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
# Function to gather optimization starting information for the Boston Harbor demographic model
def instantiate_boston_harbor():
	# Set upper and lower bounds for the list of parameter choices (in order nu1, nu2, T)
	upper_bound = [50,50,10]
	lower_bound = [1e-4,1e-4,0]
	# Initial guess at parameter values (nu1, nu2, T)
	p0 = [1,1,1]
	# Return staring information
	return upper_bound,lower_bound,p0
# Function to compute parameter values for the Boston Harbor demographic model
def params_boston_harbor(theta, mu, seqL, popt):
	# Use the sequence length, along with mu, to calculate the ancestral size
	nuA = float(theta)/float(4 * mu * seqL)
	# Use this ancestral size to calculate the nu1, nu2, and T parameter values 
	nu1 = float(popt[0]) * float(nuA)
	nu2 = float(popt[1]) * float(nuA)
	T = float(popt[2]) * float(nuA) * 2
	# Return these computed parameter values
	return [nuA,nu1,nu2,T,float(theta)]

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
# Function to gather optimization starting information for the split with mig demographic model
def instantiate_gulf_islands_ivi():
	# Set upper and lower bounds for the list of parameter choices (nu1, nu2, T, m)
	upper_bound = [50,50,10,20]
	lower_bound = [1e-4,1e-4,0,0]
	# Initial guess at parameter values (nu1, nu2, T, m)
	p0 = [1,1,1,1]
	# Return staring information
	return upper_bound,lower_bound,p0
# Function to compute parameter values for the split with mig demographic model
def params_gulf_islands_ivi(theta, mu, seqL, popt):
	# Use the sequence length, along with mu, to calculate the ancestral size
	nuA = float(theta)/float(4 * mu * seqL)
	# Use this ancestral size to calculate the nu1, nu2, T, m parameter values 
	nu1 = float(popt[0]) * float(nuA)
	nu2 = float(popt[1]) * float(nuA)
	T = float(popt[2]) * float(nuA) * 2
	m = float(popt[3])/(2 * float(nuA)) 
	# Return these computed parameter values
	return [nuA,nu1,nu2,T,m,float(theta)]

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
# Function to gather optimization starting information for the split_two_epoch demographic model with mig
def instantiate_gulf_islands_ivm():
	# Set upper and lower bounds for the list of parameter choices (in order nu1E1, nu2E1, nu1E2, nu2E2, TE1, TE2, mE1, mE2)
	upper_bound = [50,50,50,50,10,10,20,20]
	lower_bound = [1e-4,1e-4,1e-4,1e-4,0,0,0,0]
	# Initial guess at parameter values (nu1E1, nu2E1, nu1E2, nu2E2, TE1, TE2, mE1, mE2)
	p0 = [1,1,1,1,1,1,1,1]
	# Return staring information
	return upper_bound,lower_bound,p0
# Function to compute parameter values for the split_two_epoch with mig demographic model
def params_gulf_islands_ivm(theta, mu, seqL, popt):
	# Use the sequence length, along with mu, to calculate the ancestral size
	nuA = float(theta)/float(4 * mu * seqL)
	# Use this ancestral size to calculate the nu1E1, nu2E1, nu1E2, nu2E2, TE1, TE2, mE1, mE2 parameter values 
	nu1E1 = float(popt[0]) * float(nuA)
	nu2E1 = float(popt[1]) * float(nuA)
	nu1E2 = float(popt[2]) * float(nuA)
	nu2E2 = float(popt[3]) * float(nuA)
	TE1 = float(popt[4]) * float(nuA) * 2
	TE2 = float(popt[5]) * float(nuA) * 2
	mE1 = float(popt[6])/(2 * float(nuA))
	mE2 = float(popt[7])/(2 * float(nuA))
	# Return these computed parameter values
	return [nuA,nu1E1,nu2E1,nu1E2,nu2E2,TE1,TE2,mE1,mE2,float(theta)]

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
	# This is the demographic model we wish to fit 
	parser.add_argument("--model")
	# This is the optimization run number (i.e., for CHTC purposes this will index the run ID)
	parser.add_argument("--opt_num")
	# This is the mutation rate that we assume to convert parameter estimates
	parser.add_argument("--mut")
	# This is the effective sequence length that we assume (i.e., the total length of sequence from which variants *could* have been called) to convert parameter estimates
	parser.add_argument("--L")
	# Indicates whether to add a parameter for ancestral state misidentification
	parser.add_argument("--p_misid", action='store_true', default=False)

	# Parse args
	args = parser.parse_args()
	return(args)

#-------------------------
# Run the main function
#-------------------------

if __name__ == '__main__':
	main()



	
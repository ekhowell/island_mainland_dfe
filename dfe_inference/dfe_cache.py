#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# OVERVIEW
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Purpose: the purpose of this script is to cache expected FS for a collection of gammas under models of shared or independent selection
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Run as: python dfe_cache.py 
# --prefix [prefix to use for output/population(s) name(s)] 
# --fs [fs file] 
# --shared_sel_model [model name for shared selection] 
# --ind_sel_model [model name for independent selection] 
# --demog_params [demographic parameters]
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#--------------------------
# Import required packages
#--------------------------
import dadi
import dadi.DFE as DFE
import random
import numpy as np
import matplotlib.pyplot as pyplot
import argparse
import dill as pickle

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
	# Create a variable to hold the model ID for shared gammas
	shared_sel_model = args.shared_sel_model
	# Create a variable to hold the model ID for independent gammas
	ind_sel_model = args.ind_sel_model
	# Create a variable to hold the fixed parameter estimates as a list of floats
	demog_params = [float(p) for p in args.demog_params]

	#------------------------
	# Parameterizing cache
	#------------------------

	# Need to determine what the sample size is based on the input fs
	ns = fs.sample_sizes
	print("Sample sizes:")
	print(ns)

	# Need to define the grid points for optimization
	pts = [x for x in range(0, len(ns)*10, 10)] + (max(ns)*2)
	print("Grid points:")
	print(pts)

	# Wrap our demographic + selection function based on user input
	print("Wrapping models of shared and independent selection + demographic history")

	# Independent selection models
	if ind_sel_model=="boston_harbor_ind_sel":
		ind_sel_model_func = boston_harbor_ind_sel
	elif ind_sel_model=="gulf_islands_ivm_ind_sel":
		ind_sel_model_func = gulf_islands_ivm_ind_sel
	elif ind_sel_model=="gulf_islands_ivi_ind_sel":
		ind_sel_model_func = gulf_islands_ivi_ind_sel	

	# Shared selection models
	if shared_sel_model=="boston_harbor_shared_sel":
		shared_sel_model_func = boston_harbor_shared_sel
	elif shared_sel_model=="gulf_islands_ivm_shared_sel":
		shared_sel_model_func = gulf_islands_ivm_shared_sel
	elif shared_sel_model=="gulf_islands_ivi_shared_sel":
		shared_sel_model_func = gulf_islands_ivi_shared_sel

	# Wrap the demographic model in a function that utilizes grid points which increases dadi's ability to more accurately generate a model frequency spectrum.
	shared_sel_model_ex = dadi.Numerics.make_extrap_func(shared_sel_model_func)
	ind_sel_model_ex = dadi.Numerics.make_extrap_func(ind_sel_model_func)

	#----------------------------
	# Generating 1D and 2D cache
	#----------------------------

	print("Generating 1D cache (for shared selection)...")
	# Generating 1D cache
	s1 = DFE.Cache1D(demog_params, ns, shared_sel_model_ex, pts=pts, gamma_pts=100, gamma_bounds=[1e-4, 2000], verbose=False, additional_gammas=[0, 1, 2, 3, 4, 5], cpus=20)
	# Pickle the 1D (shared selection) cache
	fid = open(prefix + '_' + shared_sel_model + '_1d_cache.bpkl', 'wb')
	pickle.dump(s1, fid, protocol=2)
	fid.close()

	print("Generating 2D cache (for independent selection)...")
	# Generating 2D cache
	s2 = DFE.Cache2D(demog_params, ns, ind_sel_model_ex, pts=pts, gamma_pts=100, gamma_bounds=[1e-4, 2000], verbose=False, additional_gammas=[0, 1, 2, 3, 4, 5], cpus=20)
	# Pickle the 2D (shared selection) cache
	fid = open(prefix + '_' + ind_sel_model + '_2d_cache.bpkl', 'wb')
	pickle.dump(s2, fid, protocol=2)
	fid.close()

#----------------------------------
# Specifying the selection model
#----------------------------------

# Boston Harbor demographic model with shared selection
def boston_harbor_shared_sel(params, ns, pts):
	"""
	Parameter values:
	1. nu1 = the ratio of population 1 size to the ancestral population size
	2. nu2 = the ratio of population 2 size to the ancestral population size
	3. T = the time since the population split
	4. gamma = selection coefficient shared across the ancestral and contemporary populations
	Overview: the ancestral population splits into population 1 with size nu1 and population 2 with size nu2 T generations from the present.
	"""
	# Initialize the parameter values
	nu1, nu2, T, gamma = params	
	# Specify the density of the grid over which phi is evaluated
	xx = dadi.Numerics.default_grid(pts)
	# Ancestral dynamics #
	# Create a initial 1D phi that represents our equilibrium ancestral population with a relative population size of 1 and sel coefficient of gamma
	phi = dadi.PhiManip.phi_1D(xx, gamma = gamma)
	# Epoch 1 dynamics #
	# Split the ancestral population into two nascent populations
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	# After the splitting, the populations take constant sizes nu1 and nu2 for T generations
	phi = dadi.Integration.two_pops(phi, xx, T, nu1, nu2, gamma1 = gamma, gamma2 = gamma)
	# Return the frequency spectrum obtained from the model
	fs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))
	return fs

# Boston Harbor demographic model with indepedent selection
def boston_harbor_ind_sel(params, ns, pts):
	"""
	Parameter values:
	1. nu1 = the ratio of population 1 size to the ancestral population size
	2. nu2 = the ratio of population 2 size to the ancestral population size
	3. T = the time since the population split
	4. gamma_m = selection coefficient shared across the ancestral population and contemporary population 2 (World's End)
	5. gamma_i = selection coefficient in contemporary population 1 (Bumpkin or Peddocks)
	Overview: the ancestral population splits into population 1 with size nu1 and population 2 with size nu2 T generations from the present.
	"""
	# Initialize the parameter values
	nu1, nu2, T, gamma_m, gamma_i = params	
	# Specify the density of the grid over which phi is evaluated
	xx = dadi.Numerics.default_grid(pts)
	# Ancestral dynamics #
	# Create a initial 1D phi that represents our equilibrium ancestral population with a relative population size of 1 and sel coefficient of gamma
	phi = dadi.PhiManip.phi_1D(xx, gamma = gamma_m)
	# Epoch 1 dynamics #
	# Split the ancestral population into two nascent populations
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	# After the splitting, the populations take constant sizes nu1 and nu2 for T generations
	phi = dadi.Integration.two_pops(phi, xx, T, nu1, nu2, gamma1 = gamma_i, gamma2 = gamma_m)
	# Return the frequency spectrum obtained from the model
	fs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))
	return fs

# Define the selection model for island-mainland 2D comparisons in the Gulf Islands with shared selection between island and mainland
def gulf_islands_ivm_shared_sel(params, ns, pts):
	"""
	Parameter values:
	1. nu1E1 = the ratio of population 1 size to the ancestral population size during first epoch (mainland ancestor)
	2. nu2E1 = the ratio of population 2 size to the ancestral population size during the first epoch (island ancestor)
	3. nu1E2 = the ratio of population 1 size to the ancestral population size during second epoch (contemporary mainland)
	4. nu2E2 = the ratio of population 2 size to the ancestral population size during the second epoch (contemporary island)
	5. TE1 = the duration of the first epoch
	6. TE2 = the duration of the second epoch
	7. mE1 = the migration rate during the first epoch
	8. mE2 = the migration rate during the second epoch
	9. gamma = selection coefficient shared across all populations
	Overview: the ancestral population splits into population 1 with size nu1E1 and population 2 with size nu2E1 and migration rate mE1. After TE1 generations, the populations enter a second epoch with sizes nu1E1 and nu2E2 that lasts for TE2 generations with migration rate mE2.
	[NOTE] all parameters are given relative to the reference size. The ancestral reference size must be calculated from theta
	"""
	# Initialize the parameter values
	nu1E1, nu2E1, nu1E2, nu2E2, TE1, TE2, mE1, mE2, gamma = params
	# Specify the density of the grid over which phi is evaluated
	xx = dadi.Numerics.default_grid(pts)
	# Ancestral dynamics #
	# Create a initial 1D phi that represents our equilibrium ancestral population with a relative population size of 1
	phi = dadi.PhiManip.phi_1D(xx, gamma=gamma)
	# Epoch 1 dynamics #
	# Split the ancestral population into two nascent populations that represent the mainland ancestor and the island ancestor
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	# After the splitting, the populations take constant sizes nu1E1 and nu2E1 for TE1 generations with symmetric migration rate mE1
	phi = dadi.Integration.two_pops(phi, xx, TE1, nu1E1, nu2E1, m12=mE1, m21=mE1, gamma1=gamma, gamma2=gamma)
	# Epoch 2 dynamics #
	# The populations then take constant sizes nu1E2 and nu2E2 for TE2 generations with symmetric migration rate mE2
	phi = dadi.Integration.two_pops(phi, xx, TE2, nu1E2, nu2E2, m12=mE2, m21=mE2, gamma1=gamma, gamma2=gamma)
	# Return the frequency spectrum obtained from the model
	fs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))
	return fs

# Define the demographic function for a population split followed by constant sizes for both population 1 and 2 with symmetric migration
# This corresponds to the Saturna–Pender comparison
# This function models shared selection between island populations 
def gulf_islands_ivi_shared_sel(params, ns, pts):
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
	nu1, nu2, T, m, gamma = params
	# Specify the density of the grid over which phi is evaluated
	xx = dadi.Numerics.default_grid(pts)
	# Ancestral dynamics #
	# Create a initial 1D phi that represents our equilibrium ancestral population with a relative population size of 1
	phi = dadi.PhiManip.phi_1D(xx, gamma=gamma)
	# Epoch 1 dynamics #
	# Split the ancestral population into two nascent populations 
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	# After the splitting, the populations take constant sizes nu1 and nu2 for T generations with symmetric migration rate m
	phi = dadi.Integration.two_pops(phi, xx, T, nu1, nu2, m12=m, m21=m, gamma1=gamma, gamma2=gamma)
	# Return the frequency spectrum obtained from the model
	fs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))
	return fs

# Define the selection model for island-mainland 2D comparisons in the Gulf Islands with independent selection between island and mainland that diverges in epoch 1
def gulf_islands_ivm_ind_sel(params, ns, pts):
	"""
	Parameter values:
	1. nu1E1 = the ratio of population 1 size to the ancestral population size during first epoch (mainland ancestor)
	2. nu2E1 = the ratio of population 2 size to the ancestral population size during the first epoch (island ancestor)
	3. nu1E2 = the ratio of population 1 size to the ancestral population size during second epoch (contemporary mainland)
	4. nu2E2 = the ratio of population 2 size to the ancestral population size during the second epoch (contemporary island)
	5. TE1 = the duration of the first epoch
	6. TE2 = the duration of the second epoch
	7. mE1 = the migration rate during the first epoch
	8. mE2 = the migration rate during the second epoch
	9. gamma_m = selection coefficient shared across the ancestral population, ancestral mainland population, and contemporary Maple Ridge population
	10. gamma_i = selection coefficient in the ancestral island population and contemporary island population
	Overview: the ancestral population splits into population 1 with size nu1E1 and population 2 with size nu2E1 and migration rate mE1. After TE1 generations, the populations enter a second epoch with sizes nu1E1 and nu2E2 that lasts for TE2 generations with migration rate mE2.
	[NOTE] all parameters are given relative to the reference size. The ancestral reference size must be calculated from theta
	"""
	# Initialize the parameter values
	nu1E1, nu2E1, nu1E2, nu2E2, TE1, TE2, mE1, mE2, gamma_m, gamma_i = params
	# Specify the density of the grid over which phi is evaluated
	xx = dadi.Numerics.default_grid(pts)
	# Ancestral dynamics #
	# Create a initial 1D phi that represents our equilibrium ancestral population with a relative population size of 1
	phi = dadi.PhiManip.phi_1D(xx, gamma=gamma_m)
	# Epoch 1 dynamics #
	# Split the ancestral population into two nascent populations that represent the mainland ancestor and the island ancestor
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	# After the splitting, the populations take constant sizes nu1E1 and nu2E1 for TE1 generations with symmetric migration rate mE1
	phi = dadi.Integration.two_pops(phi, xx, TE1, nu1E1, nu2E1, m12=mE1, m21=mE1, gamma1=gamma_m, gamma2=gamma_i)
	# Epoch 2 dynamics #
	# The populations then take constant sizes nu1E2 and nu2E2 for TE2 generations with symmetric migration rate mE2
	phi = dadi.Integration.two_pops(phi, xx, TE2, nu1E2, nu2E2, m12=mE2, m21=mE2, gamma1=gamma_m, gamma2=gamma_i)
	# Return the frequency spectrum obtained from the model
	fs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))
	return fs

# Define the demographic function for a population split followed by constant sizes for both population 1 and 2 with symmetric migration
# This corresponds to the Saturna–Pender comparison
# This function models independent selection between island populations 
def gulf_islands_ivi_ind_sel(params, ns, pts):
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
	nu1, nu2, T, m, gamma1, gamma2 = params
	# Specify the density of the grid over which phi is evaluated
	xx = dadi.Numerics.default_grid(pts)
	# Ancestral dynamics #
	# Create a initial 1D phi that represents our equilibrium ancestral population with a relative population size of 1
	phi = dadi.PhiManip.phi_1D(xx, gamma=gamma1)
	# Epoch 1 dynamics #
	# Split the ancestral population into two nascent populations 
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	# After the splitting, the populations take constant sizes nu1 and nu2 for T generations with symmetric migration rate m
	phi = dadi.Integration.two_pops(phi, xx, T, nu1, nu2, m12=m, m21=m, gamma1=gamma1, gamma2=gamma2)
	# Return the frequency spectrum obtained from the model
	fs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))
	return fs

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
	# This is the 1D model we wish to fit 
	parser.add_argument("--shared_sel_model")
	# This is the 2D model we wish to fit 
	parser.add_argument("--ind_sel_model")
	# These are the best-fit demographic parameter estimates (excluding theta and the ll) for the given model
	parser.add_argument("--demog_params", default=[], nargs="+")

	# Parse args
	args = parser.parse_args()
	return(args)

#-------------------------
# Run the main function
#-------------------------

if __name__ == '__main__':
	main()



	
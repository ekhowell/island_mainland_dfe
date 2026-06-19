#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# OVERVIEW
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Purpose: the purpose of this script is to generate a dadi FS object from an input VCF file
# Output: a dadi FS object
# Special inputs: --popfile is the corresponding sample manifest in the format [individual] [pop]
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Run as: 
# python create_fs.py --vcf [input VCF file] --fspre [prefix for FS file] --popfile [corresponding sample manifest] --dimensions [2 or 3]
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#---------------------------
# IMPORT PACKAGES
#---------------------------
import pandas as pd
from random import sample
import dadi
import dill as pickle 
import random
import numpy as np
import matplotlib.pyplot as pyplot
import argparse
import pylab

#----------------------------------
# DEFINE THE MAIN SCRIPT FUNCTION
#----------------------------------

def main():

	#-----------------
	# Set variables 
	#-----------------

	print("Reading in the data and intializing variables...")

	# Parse command-line arguments
	args = parse_args()

	# Variable for holding the name of the input VCF file
	vcf = args.vcf
	# This is the name of the output FS object
	fs_name = args.fspre + ".fs"
	# Variable for holding input popfile
	popfile = args.popfile
	# Dimensions
	dim = args.dim

	print("Done.")

	#-----------------------------------------
	# Constructing frequency spectrum objects
	#-----------------------------------------

	# Run function to create FS from the input VCF
	print("Creating dadi frequency spectrum object...")
	create_fs(vcf, fs_name, popfile, dim)
	print("Done.")

#---------------------
# DEFINE FUNCTIONS
#---------------------

def parse_args():
	"""
	Parse the command-line arguments
	"""
	# Call to argparse
	parser = argparse.ArgumentParser()

	# Command-line variables

	# Determines which VCF to use
	parser.add_argument('--vcf', required=True)
	# Prefix of FS files
	parser.add_argument('--fspre', required=True)
	# Determines which popfile to use
	parser.add_argument('--popfile', required=True)
	# Specifies dimensions represented in the popfile
	parser.add_argument('--dim', required=True)

	# Parse and return arguments
	args = parser.parse_args()
	return args

def create_fs(vcf, fs_name, popfile, dim):
	"""
	This function takes in a VCF dataset and creates a dadi fs object, outputting this to a loadable .fs file and plotting in a PDF

	vcf: this is the input data
	fs_name: this is the name given to the FS fileset
	popfile: this is the corresponding sample manifest in the format [individual] [pop]

	return: NULL
	"""
	# Read in the input vcf and construct dadi dictionary object
	dd = dadi.Misc.make_data_dict_vcf(vcf, popfile)

	# Read in the provided popfile as a pandas df
	popf_df = pd.read_csv(popfile, sep=' ', header=None)

	# Identify all unique populations in the popfile
	pop_ls = np.unique(popf_df.iloc[:,1].tolist())
	print("Populations detected:")
	print(pop_ls)

	# Construct a list for the samples that correspond to each unique population
	samp_ls = [popf_df[popf_df[1] == p][0].tolist() for p in pop_ls]
	print("Corresponding samples:")
	print(samp_ls)

	# Construct a list for the sample sizes to use in the projection argument
	proj_ls = [2*len(s) for s in samp_ls]
	print("Corresponding sample sizes:")
	print(proj_ls)

	# Create a frequency spectrum for the number of populations specified in the popfile
	# Don't project down samples (i.e. not considering SNPs with missing data)
	# Polarize SNPs based on the AA field in the VCF file
	fs = dadi.Spectrum.from_data_dict(dd, pop_ls, projections = proj_ls, polarized = True)

	# Write the frequency spectrum to a file
	fs.to_file(fs_name)
	fs = dadi.Spectrum.from_file(fs_name)

	# Determine the number of dimensions to plot
	if (dim==2):
		# Open the figure
		pylab.figure(figsize=(7,6))
		# Plot the folded SFS and save as a PDF
		dadi.Plotting.plot_single_2d_sfs(fs, vmin=max(fs.min(), 1e-3))
		# Save the figure
		pyplot.show()
		pyplot.savefig((fs_name + ".pdf"))
	elif(dim==3):
		# Open the figure
		pylab.figure(figsize=(7,6))
		# Plot the folded SFS and save as a PDF
		dadi.Plotting.plot_3d_spectrum(fs, vmin=max(fs.min(), 1e-3)) 
		# Save the figure
		pyplot.show()
		pyplot.savefig((fs_name + ".pdf"))

#-------------------------
# RUN THE MAIN FUNCTION
#-------------------------

if __name__ == '__main__':
	main()





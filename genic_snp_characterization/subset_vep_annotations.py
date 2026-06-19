#!/usr/bin/python
#-------------------------------------------------------------------
# The purpose of this script is to subset a collection of VEP annotations from an inference-ready (i.e., filtered, masked, etc.) VCF file
#-------------------------------------------------------------------
# Run as: subset_vcf_coding_annotations.py 
# --vcf [inference-ready vcf] 
# --vep [vep output] 
# --popfile [population assignments, tab-delimited file formatted as (sample ID) (population)]
# --category [broad-level category]
# --annotations [VEP annotations corresponding to broad-level category]
# --out [output prefix]
#-------------------------------------------------------------------
# Outputs: this script outputs a subset version of the input VCF that only contains sites with the specified annotations. It also outputs tables with allele counts computed within each population cohort.
#-------------------------------------------------------------------

#-------------------------------------------------------------------
# LIBRARIES
#-------------------------------------------------------------------

import pandas as pd 
import argparse
import random as rd
import subprocess
import warnings
import numpy as np
warnings.simplefilter(action='ignore')

#-------------------------------------------------------------------
# MAIN
#-------------------------------------------------------------------

def main():
	"""
	Run main function
	"""
	# Parse command-line arguments
	args = parse_args()

	# Assign to variables
	vcf = args.vcf
	vep = args.vep
	popfile = pd.read_csv(args.popfile, delimiter='\t', header=None, names=['sample', 'pop'])
	category = args.category
	annotations = args.annotations
	out = args.out

	# Generate a random number
	random = rd.random()
	# Create temporary file name
	temp = 'tmp_' + str(random)

	# Create the GATK-style interval list (as a temp file) and store into memory the df with filtered annotations
	annotation_df = get_annotations(vep, annotations, temp)

	# Store into memory the df with the frequencies of annotated sites
	frequencies_df = run_gatk(vcf, temp, popfile, category, annotations, out)

#-------------------------------------------------------------------
# FUNCTIONS
#-------------------------------------------------------------------

def get_annotations(vep, annotations, temp):
	"""
	Extract sites from the VEP file corresponding to the specified chromosome and annotations
	"""
	# Specify args for the subprocess command
	args = 'cat ' + vep + ' | grep -v "#" | cut -f2,7 | awk -F":" \'{print $1"\\t"$2}\' > ' + temp
	# Use the arguments above to run the subprocess command
	subprocess.run(args, shell=True, check=True, text=True)	

	print('Storing annotations')
	# Read in the temp file as a pandas dataframe and 
	df = pd.read_csv(temp, delimiter='\t', header=None, names=['chrom', 'pos', 'annotation'])
	
	print('Removing temporary files')
	# Remove temp file
	args = 'rm ' + temp
	subprocess.run(args, shell=True, check=True, text=True)

	# Filter df to only include sites that contain at least one of the annotations
	print('Filtering annotations')
	# Use join to create a conditional statement that keeps any record containing one or more of the specified annotations
	df = df.loc[df['annotation'].str.contains(('|'.join(annotations)))]
	
	# Make a GATK-style interval list
	print('Generating GATK-style interval list')
	# Get the relevant columns and compute the 'end' position
	gatk = df[['chrom', 'pos']]
	gatk['end'] = gatk['pos'] + 1
	# Paste them together
	gatk['format'] = gatk['chrom'] + ':' + gatk['pos'].astype(str) + '-' + gatk['end'].astype(str)
	# Write the output CSV to be read directly by GATK
	gatk['format'].to_csv((temp + '_intervals.list'), index=False, header=False)

	# Return the full df with annotations for later use
	return(df)

def run_gatk(vcf, temp, popfile, category, annotations, out):
	"""
	Run the GATK SelectVariants to subset the multi-population VCF by annotation and run VariantsToTable commands to get the population-specific (and combined) allele frequencies for the same subset of sites
	"""
	# Create output prefix
	output_prefix = out + '_' + category

	print('Running GATK SelectVariants to subset full multi-population VCF')
	# Run the first SelectVariants command (once for the combined populations)
	args = '~/packages/gatk-4.2.0.0/gatk SelectVariants -V ' + vcf + ' --intervals ' + temp + '_intervals.list' + ' -O ' + output_prefix + '.vcf.gz'
	print(args)
	# Use the arguments above to run the subprocess command
	subprocess.run(args, shell=True, check=True, text=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
	# subprocess.run(args, shell=True, check=True, text=True)

	print('Removing temporary files')
	# Remove temp file
	args = 'rm ' + temp + '_intervals.list'
	print(args)
	subprocess.run(args, shell=True, check=True, text=True)

	# Obtain the list of populations present in the popfile
	pop_ls = np.unique(popfile['pop'])

	# Loop through populations in popfile and output a VariantsToTable for each
	for pop in pop_ls:

		# Subset popfile samples by population, output as temp arguments file
		samples = popfile.loc[popfile['pop'] == pop, 'sample']
		samples.to_csv((temp + '_' + pop + '.args'), index=False, header=False)

		# Use this arguments file to generate a temporary VCF with SelectVariants
		print('Running GATK SelectVariants to subset VCF by ' + pop + ' individuals')

		# Run the first SelectVariants command (once for the combined populations)
		args = '~/packages/gatk-4.2.0.0/gatk SelectVariants -V ' + output_prefix + '.vcf.gz' + ' --sample-name ' + temp + '_' + pop + '.args' + ' -O ' + temp + '_' + pop + '.vcf.gz'
		print(args)
		# Use the arguments above to run the subprocess command
		subprocess.run(args, shell=True, check=True, text=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		# subprocess.run(args, shell=True, check=True, text=True)	

		# Use this temporary VCF to output a table with VariantsToTable
		print('Running GATK VariantsToTable to obtain allele frequencies for ' + pop + ' cohort')
		args = '~/packages/gatk-4.2.0.0/gatk VariantsToTable -V ' + temp + '_' + pop + '.vcf.gz' + ' -F CHROM -F POS -F REF -F ALT -F AC -F AN -F AA -O ' + pop + '_' + output_prefix + '.table'
		print(args)
		# Use the arguments above to run the subprocess command
		subprocess.run(args, shell=True, check=True, text=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		# subprocess.run(args, shell=True, check=True, text=True)	

		print('Removing temporary files')
		# Remove temp file
		args = 'rm ' + temp + '_' + pop + '.vcf.gz* ' + temp + '_' + pop + '.args'
		print(args)
		subprocess.run(args, shell=True, check=True, text=True)

def parse_args():
	"""
	Parse the command-line arguments
	"""
	# Call to argparse
	parser = argparse.ArgumentParser()

	# List arguments
	parser.add_argument('--vcf', required=True)
	parser.add_argument('--vep', required=True)
	parser.add_argument('--popfile', required=False)
	parser.add_argument('--category', required=True)
	parser.add_argument('--annotations', required=True, nargs='+')
	parser.add_argument('--out', required=True)

	# Parse and return arguments
	args = parser.parse_args()
	return args

#-------------------------------------------------------------------
# RUN
#-------------------------------------------------------------------

if __name__ == '__main__':
	main()






















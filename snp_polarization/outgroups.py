#!/usr/bin/python
#-------------------------------------------------------------------
# The purpose of this script is to output nucleotide states in an aligned set of outgroup sequences for sites that are polymorphic in the focal species
#-------------------------------------------------------------------
# Run as: python outgroups.py --vcf [input vcf] --hal [input hal] --chrom [optional] --focal [species] --outgroups [species,species,species...] --out [output prefix]
#-------------------------------------------------------------------
# Must be run with bedtools and HAL tools accessible to script
# Can activate with source ~/packages/cactus-bin-v2.9.3/venv-cactus-v2.9.3/bin/activate
#-------------------------------------------------------------------
# This returns a dataframe in the format: [chrom] [pos] [focal state] [outgroup 1 state] [outgroup 2 state] [etc] depending on the number of outgroups
# Note that the sites returned in the output have not been filtered any way (i.e., for missing data)
# If chromosome is specified, it will be prepended to the output, otherwise all SNPs will be processed
#-------------------------------------------------------------------

#-------------------------------------------------------------------
# LIBRARIES
#-------------------------------------------------------------------

import pandas as pd 
import argparse
import random as rd
import subprocess

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
	hal = args.hal
	chrom = args.chrom
	focal = args.focal
	outgroups = args.outgroups
	out = args.out

	# Generate a random number
	random = rd.random()
	# Create temporary file name
	temp = 'tmp_' + str(random)

	# Output merged VCF intervals to be passed to halSnps
	bed = make_bed(vcf, chrom, temp)

	# Capture concatenated output of halSnps
	hal_snps(bed, temp, hal, focal, outgroups, chrom, out)

#-------------------------------------------------------------------
# FUNCTIONS
#-------------------------------------------------------------------

def make_bed(vcf, chrom, temp):
	"""
	Convert the 1-based VCF into 0-based BED format, then merge nearby intervals
	"""
	# Determine whether specific chrom should be processed
	if (chrom):
		print('Processing variants for chromosome ' + chrom)
		# Specify args for the subprocess command
		args = 'zcat ' + vcf + ' | grep -v "##" | cut -f1,2 | grep "' + chrom + '" | awk \'{print $1"\\t"$2-1"\\t"$2}\' > ' + temp + '.bed'
		# args = 'zcat ' + vcf + ' | grep -v "##" | cut -f1,2 | grep "' + chrom + '" | awk \'{print $1"\\t"$2-1"\\t"$2}\' | head -n100 > ' + temp + '.bed'
	else:
		print('No chromosome specified; processing all variants in VCF')
		# Specify args for the subprocess command
		args = 'zcat ' + vcf + ' | grep -v "##" | cut -f1,2 | awk \'{print $1"\\t"$2-1"\\t"$2}\' > ' + temp + '.bed'
		# args = 'zcat ' + vcf + ' | grep -v "##" | cut -f1,2 | awk \'{print $1"\\t"$2-1"\\t"$2}\' | head -n100 > ' + temp + '.bed'
	
	print('Making BED file...')
	# Use the arguments above to run the subprocess command
	subprocess.run(args, shell=True, check=True, text=True)
	print('Done')

	print('Merging nearby records...')
	# Specify args for running bedtools via subprocess
	args = 'bedtools merge -i ' + temp + '.bed > merged_' + temp + '.bed' 
	# Run the command with subprocess
	subprocess.run(args, shell=True, check=True, text=True)
	print('Done')

	print('Storing as data frame...')
	# Store the resulting BED records into memory
	df = pd.read_csv(('merged_' + temp + '.bed'), delimiter='\t', header=None, names=['chrom', 'start', 'stop'])
	print('Removing temporary files...')
	# Remove all temporary BED files
	args = 'rm merged_' + temp + '.bed ' + temp + '.bed'
	subprocess.run(args, shell=True, check=True, text=True)
	
	return(df)

def hal_snps(bed, temp, hal, focal, outgroups, chrom, out):
	"""
	Run halSnps command
	"""
	# Create a df to store the concatenated output from halSnps
	output = pd.DataFrame()

	print('Running halSnps on BED file...')

	# Iterate through each row of the BED file
	for index, row in bed.iterrows():
		# Store the positional arguments needed by halSnps
		chrom = str(row['chrom'])
		start = str(row['start'])
		length = str(int(row['stop']) - int(row['start']))
		
		# Construct the halSnps argument for the specified BED interval
		args = 'halSnps --minSpeciesForSnp 0 --noDupes --refSequence ' + chrom + ' --start ' + start + ' --length '+ length + ' --tsv ' + temp + '.tsv ' + hal + ' ' + focal + ' ' + ','.join(outgroups)
		# Run the halSnps command, the last arg to subprocess suppresses stdout (i.e., it is neither printed nor captured)
		subprocess.run(args, shell=True, check=True, text=True, stdout=subprocess.DEVNULL)

		# Store the output as a pandas df
		df = pd.read_csv((temp + '.tsv'), delimiter='\t')

		# Write the output directly to a file rather than store in memory (decide whether to output header based on position in BED file)
		if index == 0:
			df.to_csv((chrom + '_' + out + '.snps'), index=False, sep='\t', na_rep='NaN')
		else:
			df.to_csv((chrom + '_' + out + '.snps'), index=False, sep='\t', na_rep='NaN', header=None, mode='a')
			
	print('Done')
	print('Removing temporary files...')
	# Remove temporary tsv file
	args = 'rm ' + temp + '.tsv'
	subprocess.run(args, shell=True, check=True, text=True)

	return(output)

def parse_args():
	"""
	Parse the command-line arguments
	"""
	# Call to argparse
	parser = argparse.ArgumentParser()

	# Passes VCF file to use, including path (should be gzipped)
	parser.add_argument('--vcf', required=True)
	# Passes the HAL file to use, including path
	parser.add_argument('--hal', required=True)
	# Passes the chrom (optional)
	parser.add_argument('--chrom', required=False)
	# Passes the focal species name
	parser.add_argument('--focal', required=True)
	# Passes the outgroup species names (should be a space-separated list)
	parser.add_argument('--outgroups', required=True, nargs='+')
	# Passes the prefix to use for output
	parser.add_argument('--out', required=True)

	# Parse and return arguments
	args = parser.parse_args()
	return args

#-------------------------------------------------------------------
# RUN
#-------------------------------------------------------------------

if __name__ == '__main__':
	main()



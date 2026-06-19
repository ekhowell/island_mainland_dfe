#!/usr/bin/python
#-------------------------------------------------------------------
# The purpose of this script is to format data for est-sfs analysis
#-------------------------------------------------------------------
# Run as: python format_est_sfs.py 
# --vcf [input vcf] 
# --snps [input outgroup states from outgroups.py] 
# --chrom [chrom]
# --samps [unrelated individuals in combined populations]
# --outgroups [*ordered* list of outgroups]
# --out [output prefix]
#-------------------------------------------------------------------
# Outgroup state inference must have been run on the same VCF specified in the above arguments
#-------------------------------------------------------------------
# This returns:
# 1) A .outgroups dataframe for combined pop allele counts and outgroup states
#    Formatted as [focal] [outgroup1] [outgroup2] [outgroup3]
# 2) A corresponding .sites for the dataframe above that contains info about chrom and pos (1-based VCF coordinates) of each site
#    Formatted as [chrom] [pos]
# 3) A corresponding .info file describing which outgroups and which input files were used
# Note that the allele counts in each sample will have been pruned for missing data (computed across cohorts) and for related individuals 
#-------------------------------------------------------------------
# Example: python format_est_sfs.py 
# --vcf {input.vcf.gz}
# --snps {output.snps}
# --chrom {chrom}
# --samps {unrelated.args}
# --outgroups {selected outgroups}
# --out {output}
#-------------------------------------------------------------------

#-------------------------------------------------------------------
# LIBRARIES
#-------------------------------------------------------------------

import pandas as pd 
import argparse
import random as rd
import subprocess
import warnings
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
	snps = args.snps
	chrom = args.chrom
	pops = 'combined'
	unrel = args.samps
	outgroups = args.outgroups
	out = args.out

	# Generate a random number
	random = rd.random()
	# Create temporary file name
	temp = 'tmp_' + str(random)

	# Convert the sites coordinates in the .snps file to 1-based and recode the outgroup states, write the GATK style interval list to a temp file
	outgroup_info = format_outgroups(snps, temp, outgroups)
	print(outgroup_info)

	# Run the GATK commands and extract allele count info
	print('Getting allele count information for ' + pops + ' population(s)')
	# Assign variables for given cohort
	this_pop = pops
	this_unrel = unrel

	# Extract the allele count information for the unrelated subset from the given population. Note that this will include sites that are monomorphic or missing data in the given set (which will need to be removed at the very end)
	allele_count_info = get_allele_counts(vcf, temp, this_pop, this_unrel)
	print(allele_count_info)

	print('Merging allele counts with outgroup information for ' + pops + ' population(s)')
	# Merge the dfs containing outgroup info an allele count info for the given population. Left merge is based on keys in the outgroup_info df
	df = outgroup_info.merge(allele_count_info, on=['refSequence', 'refPosition'], how='left')
	print(df)

	print('Removing sites with missing data')
	# Determine how many individuals are in the unrelated subset
	nsamps = pd.read_csv(this_unrel, header=None)
	target_an = len(nsamps) * 2
	print(target_an)

	# Filter the df to remove missing data
	df = df[df['an'] == target_an]
	print(df)

	print('Formatting allele counts for est-sfs')
	# Tabulate the number of ref and alt alleles
	df['num_alt'] = df['ac']
	df['num_ref'] = df['an'] - df['ac']
	print(df)

	# Initialize columns to hold nucleotide counts
	df['A'] = 0
	df['C'] = 0
	df['G'] = 0
	df['T'] = 0

	# Populate those columns with the corresponding ref and alt allele counts
	df.loc[(df['ref'] == 'A') | (df['ref'] == 'a'), 'A'] = df['num_ref']
	df.loc[(df['alt'] == 'A') | (df['alt'] == 'a'), 'A'] = df['num_alt']
	df.loc[(df['ref'] == 'C') | (df['ref'] == 'c'), 'C'] = df['num_ref']
	df.loc[(df['alt'] == 'C') | (df['alt'] == 'c'), 'C'] = df['num_alt']
	df.loc[(df['ref'] == 'G') | (df['ref'] == 'g'), 'G'] = df['num_ref']
	df.loc[(df['alt'] == 'G') | (df['alt'] == 'g'), 'G'] = df['num_alt']
	df.loc[(df['ref'] == 'T') | (df['ref'] == 't'), 'T'] = df['num_ref']
	df.loc[(df['alt'] == 'T') | (df['alt'] == 't'), 'T'] = df['num_alt']

	# Combine these nucleotide columns to create a single est-sfs formatted column
	df['focal'] = df['A'].astype(str) + ',' + df['C'].astype(str) + ',' + df['G'].astype(str) + ',' + df['T'].astype(str)

	print('Generating output files')

	print('Writing est-sfs formatted file to .outgroups')
	# Subset the relevant columns and write the output for the specified cohort
	outgroup_cols = ['focal'] + outgroups
	outgroup_df = df[outgroup_cols]
	outgroup_df.to_csv((chrom + '_' + this_pop + '_' + out + '.outgroups'), index=False, header=False, sep=' ')

	print('Writing corresponding site information to .sites')
	# Subset the relevant columns and write the output for the specified cohort
	site_cols = ['refSequence', 'refPosition', 'ref', 'alt', 'ac', 'an']
	site_df = df[site_cols]
	site_df.to_csv((chrom + '_' + this_pop + '_' + out + '.sites'), index=False, sep=' ')

	print('Writing input information to .info')
	info_file = open((chrom + '_' + out + '.info'), 'w')
	info_file.write('Focal VCF: ' + vcf + '\n')
	info_file.write('Chromosome: ' + chrom + '\n')
	info_file.write('Outgroups: ' + ' '.join([str(x) for x in outgroups]) + '\n')
	info_file.write('Outgroup states: ' + snps + '\n')
	info_file.close()

	print('Removing temporary files')
	# Remove remaining temp file
	args = 'rm ' + temp + '_intervals.list'
	subprocess.run(args, shell=True, check=True, text=True)	

#-------------------------------------------------------------------
# FUNCTIONS
#-------------------------------------------------------------------

def format_outgroups(snps, temp, outgroups):
	"""
	Re-format the information in the outgroup state .snps file to reflect 1-based coordinates and to conform to the est-sfs input specification
	"""
	# Read in the .snps file as a pandas df
	df = pd.read_csv(snps, delimiter='\t', header='infer')

	# Subset the relevant position and outgroup columns
	cols = ['refSequence', 'refPosition'] + outgroups
	df = df[cols]

	print('Recoding outgroups states according to est-sfs requirements')
	# Loop through each outgroup column and recode nucleotides to est-sfs's A,G,C,T specification. Missing outgroup states are coded as 0,0,0,0
	for species in outgroups:
		df.loc[(df[species]=='a') | (df[species]=='A'), species] = '1,0,0,0'
		df.loc[(df[species]=='c') | (df[species]=='C'), species] = '0,1,0,0'
		df.loc[(df[species]=='g') | (df[species]=='G'), species] = '0,0,1,0'
		df.loc[(df[species]=='t') | (df[species]=='T'), species] = '0,0,0,1'
		df.loc[df[species].isna(), species] = '0,0,0,0'

	print('Converting positions from 0-based to 1-based coordinates')
	# Also convert the refPosition column to reflect 1-based coordinates by adding 1 to the current value
	df['refPosition'] = df['refPosition'] + 1

	print('Writing sites to GATK-style interval list')
	# Get the relevant columns and compute the 'end' position
	gatk = df[['refSequence', 'refPosition']]
	gatk['refEnd'] = gatk['refPosition'] + 1
	# Paste them together
	gatk['format'] = gatk['refSequence'] + ':' + gatk['refPosition'].astype(str) + '-' + gatk['refEnd'].astype(str)
	# Write the output CSV to be read directly by GATK
	gatk['format'].to_csv((temp + '_intervals.list'), index=False, header=False)

	# Return the df with converted coordinates and recoded outgroup information
	return(df)

def get_allele_counts(vcf, temp, pop, unrel):
	"""
	Run the GATK SelectVariants and VariantsToTable commands to get the allele count information calibrated for given cohort
	"""

	print('Running GATK SelectVariants to get unrelated individuals in ' + pop)
	# Get unrelated individuals for given cohort 
	args = '~/packages/gatk-4.2.0.0/gatk SelectVariants -V ' + vcf + ' --intervals ' + temp + '_intervals.list' + ' --sample-name ' + unrel + ' -O ' + pop + '_' + temp + '.vcf.gz'
	# Use the arguments above to run the subprocess command
	subprocess.run(args, shell=True, check=True, text=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

	print('Running GATK VariantsToTable to get allele frequencies in ' + pop)
	# Get allele frequencies calibrated to the unrelated subset
	args = '~/packages/gatk-4.2.0.0/gatk VariantsToTable -V ' + pop + '_' + temp + '.vcf.gz' + ' --intervals ' + temp + '_intervals.list' + ' -F CHROM -F POS -F REF -F ALT -F AC -F AN -O ' + pop + '_' + temp + '.table'
	# Use the arguments above to run the subprocess command
	subprocess.run(args, shell=True, check=True, text=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


	# Read in the temp file as a pandas dataframe
	df = pd.read_table((pop + '_' + temp + '.table'), delimiter='\t', skip_blank_lines=True)

	# Various shenanigans to deal with the extra column pandas adds to the VariantsToTable output
	df = df.loc[:, ~df.columns.str.startswith('Unnamed')]
	df = df.set_axis(['refSequence', 'refPosition', 'ref', 'alt', 'ac', 'an'], axis=1)

	print('Removing temporary files')
	# Remove temp files
	args = 'rm ' + pop + '_' + temp + '.vcf.gz* ' + pop + '_' + temp + '.table'
	subprocess.run(args, shell=True, check=True, text=True)	

	# Return the df with the allele counts for the unrelated subset from the specified population
	return(df)

def parse_args():
	"""
	Parse the command-line arguments
	"""
	# Call to argparse
	parser = argparse.ArgumentParser()

	# List arguments
	parser.add_argument('--vcf', required=True)
	parser.add_argument('--snps', required=True)
	parser.add_argument('--chrom', required=True)
	parser.add_argument('--samps', required=True)
	parser.add_argument('--outgroups', required=True, nargs='+')
	parser.add_argument('--out', required=True)

	# Parse and return arguments
	args = parser.parse_args()
	return args

#-------------------------------------------------------------------
# RUN
#-------------------------------------------------------------------

if __name__ == '__main__':
	main()






















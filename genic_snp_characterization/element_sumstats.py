#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# OVERVIEW
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Purpose: compute summary statistics in a per- and between- population manner for each interval in an input BED file
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Run as: 
# python element_sumstats.py --vcf [compressed VCF for observed sequence data] --bed [BED-style file] --popfile [population manifest] --out [output prefix] 
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Notes:
# BED file should be formatted as:
# [chrom] [start] [stop]
# For each element
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#---------------------------
# IMPORT PACKAGES
#---------------------------
import allel
import pandas as pd
import numpy as np
import argparse
from scipy.spatial.distance import squareform
from scipy.spatial.distance import pdist
import subprocess
import os

#----------------------------------
# DEFINE THE MAIN SCRIPT FUNCTION
#----------------------------------

def main():

	#-----------------
	# Set variables 
	#-----------------

	print("--------------------------------------------------------------")
	print("Analyzing inputs...")
	print("--------------------------------------------------------------")

	print("Intializing variables...")

	# Parse command-line arguments
	args = parse_args()
	# VCF to analyze
	invcf = args.vcf
	# BED-style file with elements
	bed = args.bed
	# Population manifest
	popfile = args.popfile
	# Output file prefix
	out = args.out

	print("Reading population manifest...")

	# Read in the popfile
	popfile = pd.read_csv(popfile, sep=" ", header=None)
	print(popfile)
	# Identify all unique populations represented in file, store names as list
	pop_ls = np.unique(popfile.iloc[:, 1].tolist())
	print("Found population:")
	print(pop_ls)
	# Identify all unique samples represented in file, store names as list
	samp_ls = np.unique(popfile.iloc[:, 0].tolist())
	print("Found samples:")
	print(samp_ls)

	print("Reading BED-style intervals...")

	# Read in the BED file
	bed = pd.read_csv(bed, sep="\t", header=None, names=['chrom', 'start', 'stop'])

	print("Creating dataframes to hold results...")

	# Create dataframes to hold results for each segment
	pi_df_ls = []
	tajD_df_ls = []
	fst_df_ls = []

	print("Iterating over element flanks...")

	# Iterate over the rows of the BED file
	for index, row in bed.iterrows():

		# Assign variables for this particular element flank
		chrom = row['chrom']
		start = row['start']
		stop = row['stop']

		# Put the flank in the correct format for the read_vcf function
		region = str(chrom) + ':' + str(start) + '-' + str(stop)
		print(region)

		# Read in the VCF file for *all* samples
		vcf = allel.read_vcf(invcf, fields='*', region=region)

		# Check that this chunk of the VCF exists (not sure why this wasn't an issue with the regional_pi.py script)
		if vcf is not None:
			# Get a list of samples in the VCF
			sample_vec = vcf['samples']
			# Extract the genotype field
			geno_vec = vcf['calldata/GT']
			# Get the vector of chromosomes (1st col in VCF)
			chrom_vec = vcf['variants/CHROM']
			# Get the vector of positions
			pos_vec = vcf['variants/POS']
			# Get the unique chromosome represented in this subset file
			chrom_ls = np.unique(chrom_vec)

			# Create a nested dictionary that will hold the population name as outermost key, then the sample names and corresponding indices in sample_vec as nested inner keys
			pop_dict = {}

			# Loop through all populations represented in popfile to create dictionary
			for name in pop_ls:
				# Add an outer key for each population
				pop_dict[name] = {}
				# Within each population key, create a key for sample names, then populate accordingly
				pop_dict[name]["samples"] = popfile[popfile.iloc[:, 1]==name].iloc[:, 0].tolist()
				# Also create a key for sample indices, then populate accordingly using their positions in the sample_vec
				pop_dict[name]["indices"] = [sample_vec.tolist().index(indv) for indv in pop_dict[name]["samples"]]
				# Also create a key to store the corresponding subset of the geno_vec
				pop_dict[name]["geno_vec"] = geno_vec[:, pop_dict[name]["indices"]]

			# Store results in a df
			pi_df = pi(pop_ls, pop_dict, geno_vec, chrom_vec, pos_vec, start, stop)
			pi_df.fillna('NA', inplace=True)
			##print(pi_df)
			# Add this df to the list
			pi_df_ls.append(pi_df)

			# Store results in a df
			tajD_df = tajD(pop_ls, pop_dict, geno_vec, chrom_vec, pos_vec, start, stop)
			tajD_df.fillna('NA', inplace=True)
			##print(tajD_df)
			# Add this df to the list
			tajD_df_ls.append(tajD_df)

			# Store results in a df
			fst_df = fst(pop_ls, pop_dict, geno_vec, chrom_vec, pos_vec, start, stop)
			fst_df.fillna('NA', inplace=True)
			##print(fst_df)
			# Add this df to the list
			fst_df_ls.append(fst_df)

	print('Combining the results...')
	# Concatenate the list of dfs together
	pi_result = pd.concat(pi_df_ls, ignore_index=True)
	tajD_result = pd.concat(tajD_df_ls, ignore_index=True)
	fst_result = pd.concat(fst_df_ls, ignore_index=True)
	
	print('Writing output files...')
	# Write output to tab-delimited file
	pi_result.to_csv((out + '.pi_per_element'), sep='\t', index=False)
	tajD_result.to_csv((out + '.tajD_per_element'), sep='\t', index=False)
	fst_result.to_csv((out + '.fst_per_element'), sep='\t', index=False)

#---------------------
# DEFINE FUNCTIONS
#---------------------

def parse_args():
	"""
	Parse the command-line arguments
	"""
	# Call to argparse
	parser = argparse.ArgumentParser()

	# File inputs

	# Passes VCF file to use (should be gzipped)
	parser.add_argument('--vcf', required=True)
	# Passes the popfile to use (identical to dadi popfile format)
	parser.add_argument('--popfile', required=False)
	# Passes the BED file to use 
	parser.add_argument('--bed', required=False)

	# Output files

	# Passes the prefix to use for output
	parser.add_argument('--out', required=True)

	# Parse and return arguments
	args = parser.parse_args()
	return args

def make_ga(geno_vec, chrom_vec, pos_vec, this_chrom):
	"""
	Make the GenotypeArray object and the positions list for the given chromosome
	"""
	# This will store the genotype array for the given chromosome
	chrom_ga = []
	# This will store the positions vector for the given chromosome
	chrom_pv = []

	# Enumerate through genotype field, storing index (i) and genotypes (g) at the ith variant in all individuals
	for i, g in enumerate(geno_vec):
		# If this variant occurs on the chromosome of interest...
		if chrom_vec[i]==this_chrom:
			# Take the genotype of *all* individuals at this site
			chrom_ga.append(g)
			# Take the position for the variant
			chrom_pv.append(pos_vec[i])

	# Turn the collection of genotypes on the given chromosome into a GA object
	chrom_ga = allel.GenotypeArray(chrom_ga)
	return chrom_ga,chrom_pv

def tajD(pop_ls, pop_dict, geno_vec, chrom_vec, pos_vec, start, stop):
	""" 
	Computes Tajima's D for a given interval
	"""
	# Empty lists to store dfs for each population
	df_ls = []

	# Iterate through all populations specified in the popfile
	for i in range(0, len(pop_ls)):

		# Assign  population name to a variable
		pop = pop_ls[i]

		# Name of current chromosome
		this_chrom = chrom_vec[0]

		# Make the GenotypeArray object and the positions list for the given chromosome
		chrom_ga,chrom_pv = make_ga(pop_dict[pop]["geno_vec"], chrom_vec, pos_vec, this_chrom)

		# Make an allele counts array from the GenotypeArray object
		chrom_ac = chrom_ga.count_alleles()
		
		# Compute windowed tajimas d across the given interval
		tajD = allel.tajima_d(chrom_ac, start=start, stop=stop, min_sites=3)

		# Create pandas df to hold the tajD results
		tajD_df = pd.DataFrame({'chrom': [this_chrom], 'start': [start], 'stop': [stop], 'tajD': [tajD], 'population': [pop]})

		# Append this populaiton df
		df_ls.append(tajD_df)

	# Concatenate the dfs
	result = pd.concat(df_ls, ignore_index=True)

	return result

def pi(pop_ls, pop_dict, geno_vec, chrom_vec, pos_vec, start, stop):
	""" 
	Computes nucleotide diversity for a given interval
	"""
	# Empty lists to store dfs for each population
	df_ls = []

	# Iterate through all populations specified in the popfile
	for i in range(0, len(pop_ls)):

		# Assign  population name to a variable
		pop = pop_ls[i]

		# Name of current chromosome
		this_chrom = chrom_vec[0]

		# Make the GenotypeArray object and the positions list for the given chromosome
		chrom_ga,chrom_pv = make_ga(pop_dict[pop]["geno_vec"], chrom_vec, pos_vec, this_chrom)

		# Make an allele counts array from the GenotypeArray object
		chrom_ac = chrom_ga.count_alleles()
		
		# Compute windowed nucleotide diversity across the given interval
		pi = allel.sequence_diversity(chrom_pv, chrom_ac, start=start, stop=stop, is_accessible=None)

		# Create pandas df to hold the pi results
		pi_df = pd.DataFrame({'chrom': [this_chrom], 'start': [start], 'stop': [stop], 'pi': [pi], 'population': [pop]})

		# Append this populaiton df
		df_ls.append(pi_df)

	# Concatenate the dfs
	result = pd.concat(df_ls, ignore_index=True)

	return result

def fst(pop_ls, pop_dict, geno_vec, chrom_vec, pos_vec, start, stop):
	"""
	Computes Fst between each pair of populations for a given interval
	"""
	# Empty lists to store dfs for each population
	df_ls = []

	# Iterate through all pairwise comparisons of the populations specified in the popfile
	for i in range(0, len(pop_ls)):
		for j in range(i+1, len(pop_ls)):

			# Assign each population name to a variable
			pop1 = pop_ls[i]
			pop2 = pop_ls[j]

			# Name of current chromosome
			this_chrom = chrom_vec[0]

			# Make the population 1 GenotypeArray object and the positions list for the given chromosome
			chrom_ga_pop1,chrom_pv_pop1 = make_ga(pop_dict[pop1]["geno_vec"], chrom_vec, pos_vec, this_chrom)
			# Make the population 2 GenotypeArray object and the positions list for the given chromosome
			chrom_ga_pop2,chrom_pv_pop2 = make_ga(pop_dict[pop2]["geno_vec"], chrom_vec, pos_vec, this_chrom)

			# Make an allele counts array for population 1 from the corresponding GenotypeArray object
			chrom_ac_pop1 = chrom_ga_pop1.count_alleles()
			# Make an allele counts array for population 2 from the corresponding GenotypeArray object
			chrom_ac_pop2 = chrom_ga_pop2.count_alleles()

			# Compute Hudsons Fst, creating relevant df columns
			fst, windows, counts = allel.windowed_hudson_fst(chrom_pv_pop1, chrom_ac_pop1, chrom_ac_pop2, size=(stop-start), start=start, stop=stop)

			# Concatenate the pop names for output
			pop = pop1 + '_' + pop2

			# Create pandas df to hold the fst results
			fst_df = pd.DataFrame({'chrom': [this_chrom], 'start': [start], 'stop': [stop], 'fst': fst, 'population': [pop]})

			# Append this populaiton df
			df_ls.append(fst_df)

	# Concatenate the dfs
	result = pd.concat(df_ls, ignore_index=True)

	return result

#-------------------------
# RUN THE MAIN FUNCTION
#-------------------------

if __name__ == '__main__':
	main()





#!/usr/bin/python
#-------------------------------------------------------------------
# The purpose of this script is to infer ancestral states based on the outputs (and inputs) to the est-sfs program
#-------------------------------------------------------------------
# Run as: python infer_ancestral_states.py 
# --sites [site info for combined populations] 
# --pvalues [pvalue info for combined populations]
# --out [output prefix]
#-------------------------------------------------------------------

#-------------------------------------------------------------------
# LIBRARIES
#-------------------------------------------------------------------

import pandas as pd 
import argparse
import random as rd
import numpy as np

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
	pop = 'combined'
	sites = args.sites
	pvalues = args.pvalues
	out = args.out

	print("Reading in site information for " + pop + "...")

	# Read in the sites information
	df = pd.read_csv(sites, sep=' ')
	print(df)

	print("Reading in ancestral state probabilities for " + pop + "...")

	# Append the ancestral state probabilities (i.e., probability of major allele being ancestral) from the est-sfs output
	df['ancProb'] = pd.read_csv(pvalues, sep=' ', header=None, usecols=[2])
	print(df)

	print("Conducting binomial samples for " + pop + "...")

	# For each site, perform a random sample of size 1 using a probability of success equal to the ancProb for that site. 1 denotes a success, 0 denotes a failure
	df['binomialSuccess'] = np.random.binomial(1, df['ancProb'])

	# Recode the binomialSuccess variable to read 'yes' (if major allele is ancestral) or 'no' (if major allele is derived)
	df.loc[df['binomialSuccess']==0, 'binomialSuccess'] = 'no'
	df.loc[df['binomialSuccess']==1, 'binomialSuccess'] = 'yes'
	print(df)

	print("Assigning ancestral states for " + pop + "...")

	# If the ac (alt allele count) is greater than an - ac, then the alt allele is major and the ref allele is minor
	df.loc[df['an'] - df['ac'] <= df['ac'], 'major'] = df['alt']
	df.loc[df['an'] - df['ac'] <= df['ac'], 'minor'] = df['ref']

	# If the ac (alt allele count) is less than an - ac, then the ref allele is major and the alt allele is minor
	df.loc[df['an'] - df['ac'] > df['ac'], 'major'] = df['ref']
	df.loc[df['an'] - df['ac'] > df['ac'], 'minor'] = df['alt']
	print(df)

	# If binomialSuccess is yes, set the major allele as ancestral 
	df.loc[df['binomialSuccess']=='yes', 'ancestral'] = df['major']
	# And set the minor allele as derived
	df.loc[df['binomialSuccess']=='yes', 'derived'] = df['minor']
	# So then the derivedAlleleFreq is the frequency of the minor (less frequent) allele
	df.loc[(df['binomialSuccess']=='yes') & (df['an'] - df['ac'] <= df['ac']), 'derivedAlleleCount'] = df['an'] - df['ac']
	df.loc[(df['binomialSuccess']=='yes') & (df['an'] - df['ac'] > df['ac']), 'derivedAlleleCount'] = df['ac']

	# If binomialSuccess is no, set the minor allele as ancestral 
	df.loc[df['binomialSuccess']=='no', 'ancestral'] = df['minor']
	# And set the major allele as derived
	df.loc[df['binomialSuccess']=='no', 'derived'] = df['major']
	# So then the derivedAlleleFreq is the frequency of the major (more frequent) allele
	df.loc[(df['binomialSuccess']=='no') & (df['an'] - df['ac'] <= df['ac']), 'derivedAlleleCount'] = df['ac']
	df.loc[(df['binomialSuccess']=='no') & (df['an'] - df['ac'] > df['ac']), 'derivedAlleleCount'] = df['an'] - df['ac']
	print(df)

	# Determine the sample size for this cohort
	nsamp = max(df['an'])

	# Write outputs
	print("Writing combined output file...")
	df.to_csv((out + '.ancestral_states'), index=False, sep=' ')
	print("Done.")

def parse_args():
	"""
	Parse the command-line arguments
	"""
	# Call to argparse
	parser = argparse.ArgumentParser()

	# List arguments
	parser.add_argument('--sites', required=True)
	parser.add_argument('--pvalues', required=True)
	parser.add_argument('--out', required=True)

	# Parse and return arguments
	args = parser.parse_args()
	return args

#-------------------------------------------------------------------
# RUN
#-------------------------------------------------------------------

if __name__ == '__main__':
	main()

























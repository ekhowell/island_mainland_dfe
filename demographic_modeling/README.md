# Re-fitting 2D demographic models for all population pairs of Boston Harbor *P. leucopus* and Gulf Islands *P. maniculatus*
This Markdown file describes the steps taken to re-fit 2D demographic models to unfolded allele frequencies for each population pair analyzed in this study. 

## Methodological details
To ensure our demographic models were well calibrated to capture the effects of population history on allele frequency dynamics at genic SNPs, we used ∂a∂i to re-fit 2D demographic models previously inferred from the *folded* SFS by Howell et al. (2025) for Boston Harbor *P. leucopus* and Howell et al. (2026) for Gulf Islands *P. maniculatus* to the *polarized* allele frequencies we obtained here for intergenic variants. To account for potential errors in our SNP polarization, we augmented these models with an additional parameter, pmisid, which represents the proportion of variants for which the ancestral state was misidentified. 

All scripts used in the analyses described below can be found in this directory.  

## Software
All of the software used for these analyses are described in the [packages](https://github.com/ekhowell/island_mainland_dfe/tree/main/packages) directory.

## Analysis and code
Prior to inference with ∂a∂i, the input VCFs containing the polarized intergenic allele frequencies for each multi-population callset must be turned into a ∂a∂i-compatible 2D FS object for each pair of populations (i.e., creating a joint, unfolded 2D SFS). This can be acheived by running the `create_fs.py` script for each pair of populations:
```
python create_fs.py --vcf {intergenic.vcf.gz} --fspre {pop1_pop2_intergenic} --popfile {pop1_pop2.popfile} --dimensions 2
```
Here, the `intergenic.vcf.gz` contains the polarized allele frequencies for intergenic SNPs identified in a given multi-population callset. The `--fspre` flag specifies what the output `.fs` file should be named. The `--popfile` flag specifies a population manifest (formatted as `[sample] [population]`) which contains the sample name and population affiliation for the given population pair. 

The resulting `.fs` files, which contain ∂a∂i-compatible, unfolded joint SFS for intergenic SNPs for each population pair can then be used to estimate demographic parameters. The `demog_models.py` conducts this maximum-likelihood estimation starting from a single intial draw of parameter values. This script specifies three distinct 2D demographic models based on the best-fitting 2D models identified for Boston Harbor *P. leucopus* (by Howell et al. (2025)) and Gulf Islands *P. maniculatus* (by Howell et al. (2026)):
1. A model called `boston_harbor` which is relevant for all island-mainland and island-island comparisons in Boston Harbor *P. leucopus*
2. A model called `gulf_islands_ivm` which is relevant for all island-mainland comparisons in Gulf Islands *P. maniculatus*
3. A model called `gulf_islands_ivi` which is relevant for the island-island comparison in Gulf Islands *P. maniculatus*

This model-fitting script can be run as follows:
```
python fit_model.py --prefix {output} --fs {pop1_pop2_intergenic.fs} --model {boston_harbor, gulf_islands_ivi, or gulf_islands_ivm} --opt_num {replicate_id} --mut {mutation_rate} --L {effective_sequence_length} --p_misid 
```
Here, the `--mut` flag specifies the per-generation, per-base pair mutation rate for the given species. The `--L` flag specifies the effective seuence length to assume (see paper Materials and Methods for details). The `--p_misid` flag instructs the script to include an additional parameter for ancestral state misidentification. 

This outputs a `.params` file containing the best-fit parameter estimates and their corresponding likelihood. This analysis should ideally be run many times for many intitial starting parameters in order to thoroughly explore the parameter space.

Once these optimizations have been performed for a large number of replicates and the highest likelihood parameter estimates have been identified, we can use the `demog_eval.py` python script to inspect the fit of the given model parameters to the observed joint (2D) and marginal (1D) unfolded site frequency spectra of intergenic SNPs for a given population pair. This script can be run as follows:
```
python demog_eval.py --prefix {output} --fs {pop1_pop2_intergenic.fs} --model {boston_harbor, gulf_islands_ivi, or gulf_islands_ivm} --popt {param1 param2 ... paramN} --p_misid
```
Where the `--popt` flag supplies a space-separated list of best-fit parameter values (excluding theta and the corresponding likelihood) and the `--p_misid` specifies whether this list of parameter values includes a parameter for ancestral state misidentification. This returns `.png` files illustrating the fit of the model as 1D and 2D residual plots. 





























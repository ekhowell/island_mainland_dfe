# Characterizing genic variation in *P.maniculatus* and *P. leucopus*
This Markdown file describes the steps taken to identify genic variants and inspect patterns of variation across genic elements in *P. maniculatus* and *P. leucopus*.

## Methodological details
Genetic variation analyzed in this study can be partitioned into “intergenic” SNPs– which are assumed to behave in a neutral manner– and “genic” SNPs– which we assume could be experiencing direct or linked selection. Our genic category is comprised of SNPs falling within NCBI RefSeq protein-coding genes. We further classify genic SNPs into five mutation types: those falling in 5’ UTRs, 3’ UTRS, introns, synonymous exon variants, and nonsynonymous exon variants. To construct these genic partitions, we used the Ensembl Variant Effect Predictor (VEP) program to generate variant predictions for SNPs observed in each *Peromyscus* multi-population callset. 

We additionally used the UCSC Genome Browser Table Browser tool to output BED-style records for introns, 5’ UTRs, 3’ UTRs, and coding exons in each species using the NCBI RefSeq gene annotation tracks (hub_6477253_ncbiRefSeq for *P. leucopus* and hub_6502785_ncbiRefSeq for *P. maniculatus*). We restricted records to autosomal genes and used bedtools to remove intervals of repetitive DNA (and inversions in *P. maniculatus*) that were masked from our SNP callset. Using these BED intervals, we then characterized patterns of variation across each element type based on our multi-population SNP callsets. Specifically, we computed nucleotide diversity (𝜋) and Tajima's D within each *P. leucopus* and *P. maniculatus* population, as well as FST between each island-mainland population pair across each element with scikit-allel.

All scripts used in the analyses described below can be found in this directory.  

## Software
All of the software used for these analyses are described in the [packages](https://github.com/ekhowell/island_mainland_dfe/tree/main/packages) directory.

## Analysis and code
Using the NCBI command-line toolkit (available [here](https://www-ncbi-nlm-nih-gov.ezproxy.library.wisc.edu/datasets/docs/v2/command-line-tools/download-and-install/)), we downloaded the NCBI datasets for both *P. maniculatus* and *P. leucopus* that include the GFF files needed by VEP to annotate variant effects. 

The following commands download, compress, and index the required *P. maniculatus* input files:
```
datasets download genome accession GCF_003704035.1 --include gff3,rna,cds,protein,genome,seq-report
unzip ncbi_dataset.zip
cp -r ncbi_dataset/data/GCF_003704035.1/ ./
rm -r ncbi_dataset ncbi_dataset.zip README.md md5sum.txt
mv GCF_003704035.1 maniculatus_ncbi_dataset
cd maniculatus_ncbi_dataset
grep -v "#" genomic.gff | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > genomic.gff.gz
tabix -p gff genomic.gff.gz
bgzip -c GCF_003704035.1_HU_Pman_2.1.3_genomic.fna > GCF_003704035.1_HU_Pman_2.1.3_genomic.fna.gz
samtools faidx GCF_003704035.1_HU_Pman_2.1.3_genomic.fna.gz
```

And the following commands perform the same operations for the required *P. leucopus* input files:
```
datasets download genome accession GCF_004664715.2 --include gff3,rna,cds,protein,genome,seq-report
unzip ncbi_dataset.zip
cp -r ncbi_dataset/data/GCF_004664715.2/ ./
rm -r ncbi_dataset ncbi_dataset.zip README.md md5sum.txt
mv GCF_004664715.2 leucopus_ncbi_dataset
cd leucopus_ncbi_dataset
grep -v "#" genomic.gff | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > genomic.gff.gz
tabix -p gff genomic.gff.gz
bgzip -c GCF_004664715.2_UCI_PerLeu_2.1_genomic.fna > GCF_004664715.2_UCI_PerLeu_2.1_genomic.fna.gz
samtools faidx GCF_004664715.2_UCI_PerLeu_2.1_genomic.fna.gz
```

Then, we can run VEP using these input annotation and assembly files to annotate variant effects for all SNPs identified in our multi-population callsets of *P. maniculatus* and *P. leucopus*.
```
vep -i {maniculatus.vcf.gz} --gff {path}/maniculatus_ncbi_dataset/genomic.gff.gz --fasta {path}/maniculatus_ncbi_dataset/GCF_003704035.1_HU_Pman_2.1.3_genomic.fna.gz --tab -o {maniculatus.annotations}
vep -i {leucopus.vcf.gz} --gff {path}/leucopus_ncbi_dataset/genomic.gff.gz --fasta {path}/leucopus_ncbi_dataset/GCF_004664715.2_UCI_PerLeu_2.1_genomic.fna.gz --tab -o {leucopus.annotations}
```
These commands produce large files (>40 GB) in BED format containing the annotations themselves (`.annotations`), along with HTML files that break down the distribution of annotations. 

We classified these variant effect predictions into five genic mutation types (introns, 5’ UTR, 3’ UTR, synonymous, and nonsynonymous) based on the Sequence Ontology hierarchy (Eilbeck et al. 2005; Cunningham et al. 2015). Intron variants include SNPs annotated as “intron_variant” (SO:0001627), “splice_donor_variant” (SO:0001575), or “splice_acceptor_variant” (SO:0001574). 5’ UTR variants are those annotated as “5_prime_UTR_variant” (SO:0001623). 3’ UTR variants are those annotated as “3_prime_UTR_variant” (SO:0001624). Synonymous variants include SNPs annotated as “synonymous_variant” (SO:0001819), “stop_retained_variant” (SO:0001567), or “start_retained_variant” (SO:0002019). Nonsynonymous variants include SNPs annotated as “nonsynonymous_variant” (SO:0001992), “stop_lost” (SO:0001578), “start_lost” (SO:0002012), “missense_variant” (SO:0001583), or “stop_gained” (SO:0001587). 

In order to partition input VCFs according to these classifications, we use the python script, `subset_vep_annotations.py`, to extract sites falling into the five genic mutation types above based on the VEP output. This script can be run according to the following command:
```
python subset_vep_annotations.py --vcf {input.vcf.gz} --vep {species.annotations} --popfile {samples.popfile} --category {intron_variant, 5_prime_UTR_variant, 3_prime_UTR_variant, synonymous_variant, or nonsynonymous_variant}   --annotations {corresponding parent and child annotations} --out {output}
```
Where the `--popfile` flag specifies a file containing the samples in the multi-population VCF file that should be output. This script outputs a subset version of the input VCF that only contains sites with the specified annotations. It also outputs tables with allele counts computed within each population cohort. 

We also used BED-style interval files for introns, 5’ UTRs, 3’ UTRs, coding exons, and intergenic regions in each species (described in the Materials and Methods section of the paper) to characterize patterns of variation across each element type. 

The python script `element_sumstats.py` uses scikit-allel to compute nucleotide diversity (𝜋) and Tajima's D within each population contained in the multi-population *P. leucopus* and *P. maniculatus* callsets, as well as FST between each island-mainland population pair. These computations are performed across each element specified in the input BED files for introns, 5’ UTRs, 3’ UTRs, coding exons, and intergenic regions. This script can be run with the following command:
```
python element_sumstats.py --vcf {input.vcf.gz} --popfile {input.popfile} --bed {input.bed} --out {output}
```
Here, the `input.popfile` specifies the population affiliation for the samples present in the `input.vcf.gz`, and the `input.bed` contains the intervals for the given element type. This outputs an `output.pi_per_element`, `output.tajD_per_element`, and `output.fst_per_element` dataframe that contains results for all population/population pairs at each element defined in the BED file input. 









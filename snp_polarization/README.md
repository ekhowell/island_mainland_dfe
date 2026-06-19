# Polarizing SNPs in *P.maniculatus* and *P. leucopus*
This Markdown file describes the steps taken to polarize single nucleotide polymorphisms identified in whole genome, multi-population callsets of *P. maniculatus* and *P. leucopus*.

## Methodological details
To gain increased resolution for downstream SFS-based analyses, we polarized SNPs in both our *P. leucopus* and *P. maniculatus* callsets, which enables us to distinguish between high and low-frequency derived variation. To polarize SNPs, we employed the approaches of Keightley et al. (2016) and Keightley and Jackson (2018) as implemented in the est-sfs program. To guide the selection of appropriate outgroup species and to determine outgroup nucleotide states at polymorphic sites, we conducted whole-genome alignments of all *Peromyscus* NCBI RefSeq assemblies available at the time of this analysis using progressive Cactus.

All scripts used in the analyses described below can be found in this directory.  

## Software
All of the software used for these analyses are described in the packages [TO DO] directory.

## Analysis and code
Prior to ancestral state inference with est-sfs, we conducted whole genome, reference-free alignments of current *Peromyscus* NCBI RefSeq assemblies. The following commands download and rename each reference assembly from NCBI's FTP site.
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/704/035/GCF_003704035.1_HU_Pman_2.1.3/GCF_003704035.1_HU_Pman_2.1.3_genomic.fna.gz
mv GCF_003704035.1_HU_Pman_2.1.3_genomic.fna.gz maniculatus.fna.gz

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/704/135/GCA_003704135.2_HU_Ppol_1.3.3/GCA_003704135.2_HU_Ppol_1.3.3_genomic.fna.gz
mv GCA_003704135.2_HU_Ppol_1.3.3_genomic.fna.gz polionotus.fna.gz

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/664/715/GCF_004664715.2_UCI_PerLeu_2.1/GCF_004664715.2_UCI_PerLeu_2.1_genomic.fna.gz
mv GCF_004664715.2_UCI_PerLeu_2.1_genomic.fna.gz leucopus.fna.gz

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/949/786/415/GCF_949786415.1_PerEre_H2_v1/GCF_949786415.1_PerEre_H2_v1_genomic.fna.gz
mv GCF_949786415.1_PerEre_H2_v1_genomic.fna.gz eremicus.fna.gz

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/007/827/085/GCF_007827085.1_ASM782708v3/GCF_007827085.1_ASM782708v3_genomic.fna.gz
mv GCF_007827085.1_ASM782708v3_genomic.fna.gz californicus.fna.gz

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/902/168/415/GCA_902168415.1_Pmel_10x_v1/GCA_902168415.1_Pmel_10x_v1_genomic.fna.gz
mv GCA_902168415.1_Pmel_10x_v1_genomic.fna.gz melanophrys.fna.gz

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/902/168/425/GCA_902168425.1_Patt_10x_v1/GCA_902168425.1_Patt_10x_v1_genomic.fna.gz
mv GCA_902168425.1_Patt_10x_v1_genomic.fna.gz attwateri.fna.gz

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/902/168/325/GCA_902168325.1_Pnud_10x_v1/GCA_902168325.1_Pnud_10x_v1_genomic.fna.gz
mv GCA_902168325.1_Pnud_10x_v1_genomic.fna.gz nudipes.fna.gz

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/902/168/405/GCA_902168405.1_Pazt_10x_v1/GCA_902168405.1_Pazt_10x_v1_genomic.fna.gz
mv GCA_902168405.1_Pazt_10x_v1_genomic.fna.gz aztecus.fna.gz
```

Prior to alignment with progressive Cactus, these assemblies must be soft-masked with RepeatMasker. The following commands perform soft-masking (i.e., converting masked bases to lowercase letters) of the input reference assemblies using default RepeatMasker parameters, the rmblast engine, and the Dfam database.
```
~/packages/RepeatMasker/RepeatMasker -engine rmblast -species 'Peromyscus' -xsmall maniculatus.fna.gz

~/packages/RepeatMasker/RepeatMasker -engine rmblast -species 'Peromyscus' -xsmall polionotus.fna.gz

~/packages/RepeatMasker/RepeatMasker -engine rmblast -species 'Peromyscus' -xsmall leucopus.fna.gz

~/packages/RepeatMasker/RepeatMasker -engine rmblast -species 'Peromyscus' -xsmall eremicus.fna.gz

~/packages/RepeatMasker/RepeatMasker -engine rmblast -species 'Peromyscus' -xsmall californicus.fna.gz

~/packages/RepeatMasker/RepeatMasker -engine rmblast -species 'Peromyscus' -xsmall melanophrys.fna.gz

~/packages/RepeatMasker/RepeatMasker -engine rmblast -species 'Peromyscus' -xsmall nudipes.fna.gz

~/packages/RepeatMasker/RepeatMasker -engine rmblast -species 'Peromyscus' -xsmall attwateri.fna.gz

~/packages/RepeatMasker/RepeatMasker -engine rmblast -species 'Peromyscus' -xsmall aztecus.fna.gz
```

Progressive Cactus also requires a guide tree, which specifies the relationships of input species. To construct this guide tree, we ran mashtree on the masked *Peromyscus* assemblies generated above.
```
mashtree --numcpus 12 *_masked.fna > peromyscus.dnd
```
Which generates the following Newick specification:
```
(attwateri_masked:0.01390,aztecus_masked:0.01695,((((polionotus_masked:0.00814,maniculatus_masked:0.00807):0.00486,leucopus_masked:0.01383):0.01537,(eremicus_masked:0.01620,californicus_masked:0.01608):0.00781):0.00758,(nudipes_masked:0.01378,melanophrys_masked:0.01290):0.00384):0.00170);
```

We use this guide tree, along with paths to the soft-masked *Peromyscus* reference assemblies as input to progressive Cactus:
```
(((attwateri:0.01390,aztecus:0.01695):0.00758,(nudipes:0.01378,melanophrys:0.01290):0.00384):0.00170,(((polionotus:0.00814,maniculatus:0.00807):0.00486,leucopus:0.01383):0.01537,(eremicus:0.01620,californicus:0.01608):0.00781):0.0);
maniculatus {path_to_softmasked_assembly}
polionotus {path_to_softmasked_assembly}
leucopus {path_to_softmasked_assembly}
eremicus {path_to_softmasked_assembly}
californicus {path_to_softmasked_assembly}
melanophrys {path_to_softmasked_assembly}
nudipes {path_to_softmasked_assembly}
attwateri {path_to_softmasked_assembly}
aztecus {path_to_softmasked_assembly}
```
This information is then added to a file called `peromyscus.txt`.

Progressive Cactus can then be run with the following command to generate a 9-way alignment in HAL format called `peromyscus.hal`.
```
cactus jobstore peromyscus.txt peromyscus.hal --maxCores 18 --workDir ./tempdir
```


Looks like everything is running as-expected with the following revised seqFile:
```
(((attwateri:0.01390,aztecus:0.01695):0.00758,(nudipes:0.01378,melanophrys:0.01290):0.00384):0.00170,(((polionotus:0.00814,maniculatus:0.00807):0.00486,leucopus:0.01383):0.01537,(eremicus:0.01620,californicus:0.01608):0.00781):0.0);
maniculatus /mnt/sas0/AD/ekhowell/projects/peromyscus_alignments/reference_genomes/maniculatus_masked.fna
polionotus /mnt/sas0/AD/ekhowell/projects/peromyscus_alignments/reference_genomes/polionotus_masked.fna
leucopus /mnt/sas0/AD/ekhowell/projects/peromyscus_alignments/reference_genomes/leucopus_masked.fna
eremicus /mnt/sas0/AD/ekhowell/projects/peromyscus_alignments/reference_genomes/eremicus_masked.fna
californicus /mnt/sas0/AD/ekhowell/projects/peromyscus_alignments/reference_genomes/californicus_masked.fna
melanophrys /mnt/sas0/AD/ekhowell/projects/peromyscus_alignments/reference_genomes/melanophrys_masked.fna
nudipes /mnt/sas0/AD/ekhowell/projects/peromyscus_alignments/reference_genomes/nudipes_masked.fna
attwateri /mnt/sas0/AD/ekhowell/projects/peromyscus_alignments/reference_genomes/attwateri_masked.fna
aztecus /mnt/sas0/AD/ekhowell/projects/peromyscus_alignments/reference_genomes/aztecus_masked.fna
```
The 9-way alignment is currently running on the screen session `cactus`. 

Running the basic halStats summary on the resulting alignment:
```
halStats peromyscus.hal
```
Returns the following:
```
hal v2.2
(((attwateri:0.0139,aztecus:0.01695)Anc3:0.00758,(nudipes:0.01378,melanophrys:0.0129)Anc4:0.00384)Anc1:0.0017,(((polionotus:0.00814,maniculatus:0.00807)Anc7:0.00486,leucopus:0.01383)Anc5:0.01537,(eremicus:0.0162,californicus:0.01608)Anc6:0.00781)Anc2:0)Anc0;

GenomeName, NumChildren, Length, NumSequences, NumTopSegments, NumBottomSegments
Anc0, 2, 2237686589, 1368, 0, 9440744
Anc1, 2, 2287523271, 3669, 11066391, 31669062
Anc3, 2, 2294839481, 16140, 31588038, 40895915
attwateri, 0, 2494792707, 171913, 42545500, 0
aztecus, 0, 2488400533, 202228, 42628764, 0
Anc4, 2, 2298391372, 15782, 31372773, 39678856
nudipes, 0, 2440382953, 179791, 40404304, 0
melanophrys, 0, 2403713740, 178386, 40162112, 0
Anc2, 2, 2300130048, 1027, 11451749, 40668467
Anc5, 2, 2327442250, 1321, 39701907, 48586588
Anc7, 2, 2373283217, 1350, 48309641, 37641436
polionotus, 0, 2645110392, 5292, 37183138, 0
maniculatus, 0, 2512439762, 8524, 37766386, 0
leucopus, 0, 2475180836, 1857, 50142369, 0
Anc6, 2, 2343812701, 959, 40503343, 54026061
eremicus, 0, 2916563079, 4297, 53682969, 0
californicus, 0, 2467923528, 237, 53834032, 0
```

We then used the `halSNPs` utility in the HAL toolkit to obtain outgroup states at sites that were observed to be polymorphic in our multi-population *P. maniculatus* and *P. leucopus* callsets. Specifically, we use the `--refSequence`, `--start`, and `--length` arguments in `hapSnps` which specify the chromosome, start position, and length of interval over which to output nucleotide states in the focal and outgroup species. This analysis is performed by the python wrapper script `outgroups.py`, which performs the following steps:
1. Takes as arguments an input VCF, the chromosome to analyze, the HAL alignment file, name of focal species, name of outgroups to get nucleotide states for, and the prefix to use for outputs.
2. Extracts SNP calls from the chromosome (if specified), converts them the BED format, then merges nearby intervals to create a streamlined BED.
3. For each entry in the BED file, passes the refSequence, start, and length to halSnps (along with other arguments).
4. Converts output to a pandas df, then appends it to a growing df for either the specified chromosome or the whole genome

For each multi-population callset, we run the script on a per-chromosome basis to generate the nucleotide state in each outgroup at the sites present in the multi-population VCF. Below is an example call to this script:
```
python outgroups.py --vcf {input.vcf.gz} --hal peromyscus.hal --chrom {chromosome} --focal {maniculatus or leucopus} --outgroups leucopus polionotus eremicus californicus melanophrys nudipes attwateri aztecus --out {output.snps}
```
The resulting `output.snps` file contains the following columns:
```
[chrom] [pos] [focal state] [outgroup 1 state] [outgroup 2 state] [etc] depending on the number of outgroups
```
Where positions are specificed in 0-based coordinates. 

Once these outgroups states are determined, we have to format input files according to est-sfs specifications (as described in the program documentation):
> If there are three outgroups, there are 4 space-separated columns. The first column is for the focal species, and the next three columns are for the outgroups. Each column is a comma-separated list of the counts of the four bases in the order A, C, G, T. For example, the first line in the example data file is: `20,0,0,0 0,0,0,1 0,0,0,1 0,0,0,1` At this site, all n = 20 copies sampled in the focal species are A. In the three outgroups, a single copy has been sampled, and in each case it is T. All sites must have the same number of copies sampled in the focal species and up to one copy sampled in each outgroup. If there are missing data in any outgroup, the counts for that outgroup are encoded 0,0,0,0 . Data from polymorphic and non-polymorphic sites are analysed together.

To generate the inputs from the `input.vcf.gz` and `output.snps` files, we use a python wrapper script called `format_est_sfs.py`, which performed the following steps:
1. Formats the outgroup columns of the .snps file in the A,C,G,T format required by the method. Preserve the refSequence and refPosition column, adding 1 to the refPosition. 
2. Output the coordinate-converted positions file to a GATK-style interval list, then use this to extract the REF, ALT, AC, and AN lists from the multi-population callset.
3. Formats the focal allele counts into the A,C,G,T format.
4. Merges this focal allele count information with the outgroup information, remove sites with missing data and split the resulting filtered file into a `.outgroups` and a `.sites` file where the former is formatted according to est-sfs requirements and the latter contains the corresponding site information (much like a PLINK .map file)

This script can be run with the following command:
```
python format_est_sfs.py --vcf {input.vcf.gz} --snps {output.snps} --chrom {chrom} --samps {unrelated.args} --outgroups {selected outgroups} --out {output}
```
Where the `--outgroups` flag specifies the outgroups that should be used and the `--samps` flag specifies a file containing a list of samples in the multi-population VCF callset for which allele frequencies should be computed. 

This outputs three files:
1. An `output.outgroups` dataframe for combined allele counts in the multi-population callset and outgroup states, formatted as `[focal] [outgroup1] [outgroup2] [outgroup3]`
2. A corresponding `output.sites` file for the dataframe above that contains information about chrom and pos (1-based VCF coordinates) of each site, formatted as `[chrom] [pos]`
3. A corresponding `output.info` file describing which outgroups and which input files were used

The est-sfs program can then be run with the following command.
```
est-sfs config.txt output.outgroups seedfile.txt output.sfs output.pvalues
```
Where the `config.txt` file is formatted as follows:
```
n_outgroup 3
model 1
nrandom 10
```
To specify 3 outgroups, the Kimura 2-parameter model of nucleotide substitution, and the number of random starting value runs to use for ML estimation. And the seed file, `seedfile.txt`, contains a single random number.

The resulting outputs, `output.sfs` and `output.pvalues` contain the inferred uSFS (which we ignore in this case), and the ancestral state probabilities (which we use for downstream ancestral state inference), respectively. 

Using the `output.pvalues` file, we sample ancestral state probabilities according to a binomial distribution, determine which allele is ancestral, then add this information to the earlier `output.sites` file the corresponding derived allele and derived allele frequency.

To do this, we use the `infer_ancestral_states.py` python wrapper script, which performs the following steps:
1. Grabs all columns from the `output.sites` file and adds the third column (minus header) from the `output.pvalues` file, which correspond to the probability that the major allele is ancestral.
2. Uses binomial sampling to determine whether the major allele is derived (binomialSuccess).
3. Based on this result, adds two columns (ancestral and derived) that state the corresponding alleles based on the binomialSuccess result.
4. Adds two columns (ancestralFreq and derivedFreq) that contain the ancestral and derived frequencies to the `output.sites` file. 

This `infer_ancestral_states.py` script can be run with the following commands:
```
python infer_ancestral_states.py --sites {output.sites} --pvalues {output.pvalues} --out {output}
```
Which returns an `output.ancestral_states` file formatted as follows:
```
[refSequence] [refPosition] [ref] [alt] [ac] [an] [ancProb] [binomialSuccess] [major] [minor] [ancestral] [derived] [derivedAlleleCount]
```

Finally, the ancestral state determinations predent in the `output.ancestral_states` file can be used to annotate the AA field in the original VCF file with bcftools. 

First, we must create a compressed, tab-delimited annotation file (plus index) in the format `CHROM POS REF ALT INFO/AA` using the corresponding columns of the `output.ancestral_states` file:
```
cat output.ancestral_states | grep -v "refSequence" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$11}' > AA.tab
bgzip -c AA.tab > AA.tab.gz
tabix -s1 -b2 -e2 AA.tab.gz
```
Then, we create an INFO header line to use in the AA-annotated VCF:
```
echo '##INFO=<ID=AA,Number=1,Type=Character,Description="Ancestral allele">' > hdr.txt
```
Finally, we run the bcftools command to add both the header line and the ancestral allele annotations to the VCF:
```
bcftools annotate -a AA.tab.gz -c CHROM,POS,REF,ALT,INFO/AA -h hdr.txt -Oz -o {annotated.vcf.gz} {input.vcf.gz}
```
And remove any sites for which the ancestral allele could not be identifed:
```
bcftools view -e 'INFO/AA=="."' -Oz -o {filtered.vcf.gz} {annotated.vcf.gz}
```













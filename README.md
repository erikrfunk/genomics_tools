# Mostly python scripts to perform some basic file manipulations and statistics on fasta / VCF files
--------------------------------
# vcf2fasta
Uses VCFparser to read in a vcf and output a fasta file
Two options to deal with heterozygotes using the argument
[whichsnp]
'random': Either the ref or alt allele is chosen at random
'pair': Both alleles are written to the fasta file as a pair, with the ref allele written first
Additional arguments include:
[outfile] Name for the fasta output files
[pop] Name of a population of individuals to use. If specified, a popsFile needs to be included. See next.
[popsFile] File containing sample names in first column and population name in the second.
[individuals] Still in development, but will be used to include only select individuals from either text input or a file

# genomics_vcf_summary_stats
For now this is really only useful to calculate Watterson's theta for a single, or set of populations. This will eventually include calculations for Pi and Dxy
To calculate theta, call Wattersons_theta
[snps_file] A vcf file
[pops] A file with samples in first column and population name in the second
[pop] A population to calculate theta for. Defulat is 'None' in which theta is calculated for all populations

# VCFparser
Includes class 'Variants' that parses a vcf file and returns general attributes \
Two class functions, 'genotype_convert' and 'index_pops'
Multiple format types can be passed as arguments to genotype convert depending on the format desired
Index pops takes a txt file with population data to return the field indices for the individuals in that pop
This module is imported by most all scripts in this repo that use vcf files

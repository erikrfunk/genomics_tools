# Mostly python scripts to perform some basic file manipulations and statistics on fasta / VCF files
--------------------------------
# vcf2fasta
Uses VCFparser to read in a vcf and output a fasta file
Two options to deal with heterozygotes using the argument \
[whichsnp] \
'random': Either the ref or alt allele is chosen at random \
'pair': Both alleles are written to the fasta file as a pair, with the ref allele written first \

Additional arguments include: \
[outfile] Name for the fasta output files \
[pop] Name of a population of individuals to use. If specified, a popsFile needs to be included. See next. \
[popsFile] File containing sample names in first column and population name in the second. \
[individuals] Still in development, but will be used to include only select individuals from either text input or a file

# genomics_vcf_summary_stats
For now this is really only useful to calculate Watterson's theta for a single, or set of populations. This will eventually include calculations for Pi and Dxy
To calculate theta, call Wattersons_theta \
[snps_file] A vcf file \
[pops] A file with samples in first column and population name in the second \
[pop] A population to calculate theta for. Defulat is 'None' in which theta is calculated for all populations

# VCFparser
Includes class 'Variants' that parses a vcf file and returns general attributes \
Two class functions, 'genotype_convert' and 'index_pops'
Multiple format types can be passed as arguments to genotype convert depending on the format desired
Index pops takes a txt file with population data to return the field indices for the individuals in that pop
This module is imported by most all scripts in this repo that use vcf files

# Per-variant_Hetero_Counts
Counts the proportion of heterozygotes for each variant to look for possible paralogs

update on commands forthcoming

# sample-random-frac.py
Selects a random fraction of variants from a VCF file. This operates like the "select-random-frac"
command in GATK but uses a random seed to sample different variants each time. \
[--input] Name of input vcf \
[--fraction] Fraction of variants to sample (from 0-1) \
[--output] Defaults to 'random-variants.vcf'

# avg-depth-from-DP.py
Generating a vcf table from samtools 'mpileup' and bcftools 'call' seems not to
include individual depth information in the format field of each variant. This
prevents a few of VCFtools' filter functions from working. The samtools vcf does
include the total variant depth in the info field of each variant. Using the
--summary flag, this script generates a summary of the average depth for each variant by dividing the info field by the number of samples in the vcf table.
Instead of calling the --summary argument, min and max values can be provided
to filter over or under covered variants to remove potential copy number
variants. \
[--vcf] Name of input vcf file \
[--output] Defaults to output.depth \
[--summary] Set this to generate a summary only. No argument required. No
filtering will occur, even if filters are set. \
[--mindp] If summary flag is not set, variants below this value will be filtered \
[--maxdp] If summary flag is not set, variants above this value will be filtered

#!/usr/bin/env python3

'''
Parses a VCF file to tally the total number of heterozygotes per variant
Only argument is a vcf file

imports from VCFparser'''

import sys
from VCFparser import Variants

#Parameters
vcfPath = sys.argv[1] #vcf file to parse

#Format genotypes to count
sys.stderr.write("Formatting genotypes")
vcf = Variants(vcfPath)
genotypeLines = vcf.genotype_convert(genoformat = 'numeric')

#Tally will create a dict where the key is the variant number and value is the
#proportion of individuals genotyped at that locus that are heterozygous
sys.stderr.write("Tallying heterozygotes per variant")
heteroCount = {}
variantCount = 0

for variant in genotypeLines:
    variantCount += 1
    numGenotypedSamples = 0
    numHeterozygotes = 0
    allGenotypes = variant[9:] #Cut out the position and ref info
    for genotype in allGenotypes:
        alleles = genotype.split(sep = '/')
        if alleles[0] == alleles[1] and '.' not in alleles:
            numGenotypedSamples += 1
        elif alleles[0] != alleles[1] and '.' not in alleles:
            numGenotypedSamples += 1
            numHeterozygotes += 1
    heteroCount[variantCount] = numHeterozygotes / numGenotypedSamples

with open('Per-variant_Hetero_Counts.txt', 'w') as o:
    for i in range(1, variantCount+1):
        o.write("%s\t%s\n" % (str(i), str(heteroCount[i])))

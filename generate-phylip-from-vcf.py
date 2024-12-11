#!/usr/bin/env python
'''
Given a vcf file, this script will convert variants to phylip for the
listed taxa. If no taxon file is provided, all individuals are included.
Allele for heterozygotes is chosen at random.
'''

import sys
import argparse
import random
import gzip
from cyvcf2 import VCF

def get_args():
    parser = argparse.ArgumentParser(description='Generate concatenated sequence')
    parser.add_argument(
        '--vcf',
        required = True,
        help='VCF file'
    )
    parser.add_argument(
        '--indvfile',
        required = False,
        default = None,
        help = 'A file with sample names of individuals to include in the output'
    )
    parser.add_argument(
        '--output',
        required = False,
        default = "output",
        help = "Name of output file"
    )
    parser.add_argument(
        '--derived',
        required = False,
        default = False,
        action = 'store_true',
        help = 'If added, the alternate allele will always be used for heterozygotes'
    )
    return parser.parse_args()

def get_samples(input):
    if input[-3:] == '.gz':
        with gzip.open(input, 'rt') as file:
            for line in file:
                if line[0] == '#' and line[0:2] != '##':
                    samples = line.split()[9:]
                    return samples
    else:
        with open(input, 'r') as file:
            for line in file:
                if line[0] == '#' and line[0:2] != '##':
                    samples = line.split()[9:]
                    return samples


def check_biallelic(input_vcf):
    alternates = []
    for v in VCF(input_vcf):
        alternates.append(''.join(v.ALT))
    biallelic = True
    for allele in alternates:
        if len(allele)>1:
            biallelic = False
    return biallelic

def parse_vcf(input_vcf, sample_idx,derived):
    variant_list = []
    for v in VCF(input_vcf):
        alleles = v.genotypes[sample_idx]
        geno = str(alleles[0])+str(alleles[1])
        if geno in ['00', '.0', '0.']:
            variant_list.append(''.join(v.REF))
        elif geno in ['11','.1','1.']:
            variant_list.append(''.join(v.ALT))
        else:
            if derived == False:
                variant_list.append(''.join(random.choice([v.REF,v.ALT])))
            else:
                variant_list.append(''.join(v.ALT))
    return variant_list


def main():
    arguments = get_args()

    # Get an enumerated list of samples, and if an indv list is provided,
    # subsample it
    samples = get_samples(arguments.vcf)
    sample_idxs = [(a,b) for a,b in enumerate(samples)]

    if arguments.indvfile is not None:
        subsamples = []
        with open(arguments.indvfile, 'r') as indv_file:
            for line in indv_file:
                subsamples.append(line.strip())

        idxs = [samples.index(x) for x in subsamples]
        sample_idxs = list(zip(idxs,subsamples)) # overwrites the original indx list
    # Check that all sites are biallelic by taking the length of the alternate
    biallele_bool = check_biallelic(arguments.vcf)
    if biallele_bool == False:
        print("""
        Some alternate variants have a length greater than one.
        This might suggest a site that is more than biallelic.
        Exiting...
        """)
        exit()

    # Iterate through the sample list and add variants to the phylip
    f_out = open(str(arguments.output)+'.phylip','w')
    count=0
    for individual in sample_idxs:
        id = str(individual[1])
        idx = individual[0]
        variants = parse_vcf(arguments.vcf,idx,arguments.derived)
        if count==0:
            f_out.write('\t'.join((str(len(sample_idxs)),str(len(variants)),'\n')))
            f_out.write('\t'.join((id,''.join(variants),'\n')))
        else:
            f_out.write('\t'.join((id,''.join(variants),'\n')))
        count+=1
        print(id, ' added.', str(len(sample_idxs)-count), 'remaining.')
    f_out.close()


############################

if __name__ == '__main__':
    main()

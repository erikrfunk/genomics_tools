#!/usr/bin/env python
'''
Given a vcf file, this script will convert variants to phylip for the
listed taxa. If no taxon file is provided, all individuals are included.
Allele for heterozygotes is chosen at random unless --derived flag is used.
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
    parser.add_argument(
        '--check',
        required = False,
        default = False,
        action = 'store_true',
        help = "Performs a check to make sure all alleles are biallelic"
    )
    parser.add_argument(
        '--lowmem',
        required = False,
        default = False,
        action = 'store_true',
        help = "Writes output for one individual at a time. Memory efficient at the expense of speed"
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

def parse_vcf_lowmem(input_vcf, sample_idx, derived):
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


def parse_vcf(input_vcf, sample_idx, derived):
    var=0
    variant_list = []
    variant_tmp = []
    for v in VCF(input_vcf):
        ref = ''.join(v.REF)
        alt = ''.join(v.ALT)
        var_dict_derived = {
        '00':ref, '0-1':ref, '-10':ref,
        '11':alt, '1-1':alt, '-11':alt,
        '01':alt, '10':alt,'-1-1':'N'
        }
        var_dict = {
        '00':ref, '0-1':ref, '-10':ref,
        '11':alt, '1-1':alt, '-11':alt,
        '01':random.choice([ref,alt]), '10':random.choice([ref,alt]),'-1-1':'N'
        }
        alleles = [v.genotypes[x] for x in sample_idx]
        genos = [str(x[0])+str(x[1]) for x in alleles]
        if derived == True:
            for key,value in var_dict_derived.items():
                if key in genos:
                    genos = [i.replace(key,value) for i in genos]
        else:
            for key,value in var_dict.items():
                if key in genos:
                    genos = [i.replace(key,value) for i in genos]
        if var == 0:
            variant_tmp = genos
        else:
            if var % 50000 == 0:
                print(str(var), "variants finished.")
            if var % 10000 == 0:
                variant_list.append(variant_tmp)
                variant_tmp = genos
            else:
                variant_tmp = [x+y for x,y in zip(variant_tmp,genos)]

        var+=1
    variant_list.append(variant_tmp)
    all_vars = [''.join(x) for x in zip(*variant_list)]
    return all_vars

def main():
    arguments = get_args()
    print(arguments)
    # Get an enumerated list of samples, and if an indv list is provided,
    # subsample it
    print("Retrieving samples.")
    samples = get_samples(arguments.vcf)
    sample_idxs = [(a,b) for a,b in enumerate(samples)]

    if arguments.indvfile is not None:
        subsamples = []
        with open(arguments.indvfile, 'r') as indv_file:
            for line in indv_file:
                subsamples.append(line.strip())

        idxs = [samples.index(x) for x in subsamples]
        sample_idxs = list(zip(idxs,subsamples)) # overwrites the original indx list

    print("Writing variants for", str(len(sample_idxs)), "samples.")

    # Check that all sites are biallelic by taking the length of the alternate
    if arguments.check == True:
        print("Checking variants.")
        biallele_bool = check_biallelic(arguments.vcf)
        if biallele_bool == False:
            print("""
            Some alternate variants have a length greater than one.
            This might suggest a site that is more than biallelic.
            Exiting...
            """)
            exit()

    # Iterate through the sample list and add variants to the phylip
    if arguments.lowmem == True:
        f_out = open(str(arguments.output)+'.phylip','w')
        count=0
        for individual in sample_idxs:
            id = str(individual[1])
            idx = individual[0]
            variants = parse_vcf_lowmem(arguments.vcf,idx,arguments.derived)
            if count==0:
                f_out.write('\t'.join((str(len(sample_idxs)),str(len(variants)),'\n')))
                f_out.write('\t'.join((id,''.join(variants),'\n')))
            else:
                f_out.write('\t'.join((id,''.join(variants),'\n')))
            count+=1
            print(id, ' added.', str(len(sample_idxs)-count), 'remaining.')
        f_out.close()
    else:
        idx = [x[0] for x in sample_idxs]
        id = [x[1] for x in sample_idxs]
        variants = parse_vcf(arguments.vcf,idx,arguments.derived)
        with open(str(arguments.output)+'.phylip','w') as f_out:
            f_out.write('\t'.join((str(len(sample_idxs)),str(len(variants)),'\n')))
            [f_out.write('\t'.join([x[0],x[1],'\n'])) for x in list(zip(id,variants))]

############################

if __name__ == '__main__':
    main()

#!/usr/bin/env python
'''
Given an input of reference genome, vcf file, and a command line blast fmtout 6
file, this script constructs a multi fasta file for the locus of interest.
'''

import sys
import argparse
import gzip


def get_args():
    parser = argparse.ArgumentParser(description='Extract locus using blast result')
    parser.add_argument(
        '--ref',
        required = True,
        help = 'Reference fasta'
    )
    parser.add_argument(
        '--vcf',
        required = True,
        default='VCF file compressed with gzip'
    )
    parser.add_argument(
        '--locus',
        required = True,
        default=None,
        help = 'A command line blast result, in fmtout 6 format'
    )
    parser.add_argument(
        '--indiv',
        required = False,
        default = None,
        help = 'A file with sample names of individuals to include in the fasta'
    )
    parser.add_argument(
        '--output',
        required = False,
        default = 'output',
        help = 'Prefix for output file')
    return parser.parse_args()


def trim_vcf(vcf,indivs,scaff,start,stop):
    '''
    Trim the vcf to only include the relevant individuals,
    scaffold, and sites
    '''
    new_vcf = []
    sample_indicies = {}
    with gzip.open(vcf, 'r') as vcf_in:
        for line in vcf_in:
            if line[0:2] != '##' :
                if line[0] == '#':
                    line = line.split()
                    header = line[9:]
                    if indivs is not None:
                        sys.stderr.write("Retaining only:\n{}\n".format(indivs))
                        for individual in indivs:
                            sample_indicies[individual] = \
                            header.index(individual)+9
                    else:
                        for individual in header:
                            sample_indicies[individual] = \
                            header.index(individual)+9
                else:
                    line = line.split()

                    if (line[0] == scaff and int(line[1]) >= start and
                            int(line[1]) <= stop):
                        new_vcf.append(line)

    return new_vcf, sample_indicies


def parse_trimmed_vcf(trimmed_vcf, sample_indicies):
    '''
    Parse a trimmed vcf file saved as variable from 'trim_vcf'
    Creates a dict with individuals as keys and alleles as values for all
    variant sites
    Only one allele is included per site
    If the individual posses at least one alternate allele, that allele is written
    '''
    variant_alleles = {}
    variant_sites = []
    for line in trimmed_vcf:
        variant_sites.append(line[1])
        for sample in sample_indicies:
            if sample not in variant_alleles:
                variant_alleles[sample] = []
            vcf_index = sample_indicies[sample]
            sample_variant_info = line[vcf_index].split(':')
            if (sample_variant_info[0] == '0/0' or
                sample_variant_info[0] == './.'):
                variant_alleles[sample].append(line[3])
            elif (sample_variant_info[0] == '1/0' or
                  sample_variant_info[0] == '0/1' or
                  sample_variant_info[0] == './1' or
                  sample_variant_info[0] == '1/.' or
                  sample_variant_info[0] == '1/1'):
                variant_alleles[sample].append(line[4])
            else:
                sys.stderr.write('Variant sites may not all be biallelic\n')
    return variant_alleles, variant_sites

###############################################################################

def main():
    arguments = get_args()
    sys.stderr.write("Command has been called with the following "\
                     "arguments:\n{}\n".format(arguments))


    #Read in the reference sequence as dict with scaff name as key and seq as
    #value
    sys.stderr.write("\nReading in reference sequence...\n")
    reference = {}
    with open(arguments.ref, 'r') as ref_file:
        for line in ref_file:
            newline = line.strip()
            if newline[0] == '>':
                scaffold = newline[1:]
                reference[scaffold] = ''
            else:
                reference[scaffold] += newline


    #Get the location of the locus
    scaffold = ''
    start_pos = 0
    end_pos = 0
    with open(arguments.locus, 'r') as locus_file:
        hits = [line.strip() for line in locus_file]
        best_hit = hits[0].split()
        scaffold = best_hit[1]
        poses = [int(x) for x in best_hit[8:10]]
        start_pos = min(poses)
        end_pos = max(poses)
    sys.stderr.write("Extracting positions {}-{} of scaffold {}\n"\
                    .format(start_pos,end_pos,scaffold))

    #call the vcf trim function defined above
    #First check to see if specific indivuals have been included
    if arguments.indiv is not None:
        with open(arguments.indiv, 'r') as indiv_file:
            individuals = [line.strip() for line in indiv_file.readlines()]
    else:
        individuals = None
    trimmed_vcf, indicies = trim_vcf(arguments.vcf, individuals, scaffold,
            start_pos, end_pos)
    sys.stderr.write("{} variants found within this "\
                     "locus\n".format(len(trimmed_vcf)))


    #Parse the trimmed vcf to get genotypes/allele
    alternate_alleles, variant_pos = parse_trimmed_vcf(trimmed_vcf, indicies)

    #Make a new dict with each sample as key and the ref seq as value
    #Change variant sites in reference for each individual
    full_scaffold_dict = {}
    fout = open(arguments.output, 'w')
    for id in alternate_alleles:
        temp_seq_list = [x for x in reference[scaffold]]
        for position in range(0,len(alternate_alleles[id])):
            temp_seq_list[int(variant_pos[position])-1] = \
                    alternate_alleles[id][position]
        final_seq = temp_seq_list[start_pos-1:end_pos]
        final_seq = ''.join(final_seq)
        if poses[0] > poses[1]:
            final_seq = final_seq[::-1] #Need to reverse the direction if the blast hit was reversed
            fout.write(">{}\n{}\n".format(id,final_seq))
        else:
            fout.write(">{}\n{}\n".format(id,final_seq))
    fout.close()


###############################################################################

if __name__ == '__main__':
    main()

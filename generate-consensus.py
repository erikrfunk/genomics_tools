#!/usr/bin/env python
'''
Given an input of reference genome, and a vcf file, this script will generate
consensus fastas for the listed taxa. If no taxa or file is provided, a
consensus fasta will be created for each individual.
'''

import sys
import argparse

def get_args():
    parser = argparse.ArgumentParser(description='Generate consensus sequence')
    parser.add_argument(
        '--ref',
        required = True,
        help = 'fasta file of the reference sequence'
    )
    parser.add_argument(
        '--vcf',
        required = True,
        help='VCF file'
    )
    parser.add_argument(
        '--indv',
        required = False,
        default=None,
        help = 'A comma separated list of individuals to include'
    )
    parser.add_argument(
        '--indvfile',
        required = False,
        default = None,
        help = 'A file with sample names of individuals to include in the output'
    )
    parser.add_argument(
        '--scaffold',
        required = False,
        default = None,
        help = 'A comma separated list of scaffold(s) to construct consensus for'
    )
    parser.add_argument(
        '--scaffoldfile',
        required = False,
        default = None,
        help = 'A file with a list of scaffolds to include'
    )
    parser.add_argument(
        '--format',
        required = False,
        default = 'fasta',
        help = 'output format: fasta, or phylip'
    )
    parser.add_argument(
        '--variants',
        required = False,
        default = False,
        action = 'store_true',
        help = 'If true, only variant sites will be printed. Phylip only'
    )
    parser.add_argument(
        '--output',
        required = False,
        default = "ConsensusSeqs",
        help = "Name of output file"
    )
    return parser.parse_args()


def parse_ref(ref_seq):
    reference = {}
    genome_length = 0
    with open(ref_seq, 'r') as ref_file:
        for line in ref_file:
            newline = line.strip()
            if newline[0] == '>':
                scaffold = newline[1:]
                reference[scaffold] = ''
            else:
                reference[scaffold] += newline
                genome_length += len(newline)
    return reference, genome_length


def get_samples(input):
    with open(input, 'r') as file:
        for line in file:
            if line[0] == '#' and line[0:2] != '##':
                samples = line.split()[9:]
    return samples


def parse_vcf(input_vcf, sample, scaffold_list):
    '''
    Fetch the variants and site position for an indvidual
    '''
    variants = {}
    sites = {}
    sample_index = 0
    genos = []
    number = 0
    with open(input_vcf, 'r') as vcf:
        for line in vcf:
            if line[0:2] != "##":
                if line[0] == '#':
                    sample_index = line.split().index(sample)
                else:
                    number += 1
                    info = line.split()
                    genotype = []
                    if scaffold_list is None or info[0] in scaffold_list:
                        if info[0] in sites:
                            sites[info[0]].append(info[1])
                        else:
                            sites[info[0]] = [info[1]]
                            variants[info[0]] = []
                        genotype = info[sample_index].split(':')
                        if (genotype[0] == '0/0'):
                            variants[info[0]].append(info[3])
                        elif (genotype[0] == '1/0' or
                              genotype[0] == '0/1' or
                              genotype[0] == './1' or
                              genotype[0] == '1/.' or
                              genotype[0] == '1/1'):
                            variants[info[0]].append(info[4])
                        else:
                            variants[info[0]].append('N')

    return variants, sites, number


################################################################################

def main():
    #Parse the reference and generate sample/scaffold list
    arguments = get_args()
    sys.stderr.write("\nReading in reference sequence...\n")
    reference,whole_genome_length = parse_ref(arguments.ref)

    sample_list = []
    if arguments.indv is not None:
        sys.stderr.write("Retaining only: {}\n".format(arguments.indv))
        sample_list = arguments.indv.split(',')
    elif arguments.indvfile is not None:
        with open(arguments.indvfile, 'r') as indv_file:
            for line in indv_file:
                sample_list.append(line)
        sys.stderr.write("Retaining {} individuals\n".format(len(sample_list)))
    else:
        sample_list = get_samples(arguments.vcf)
        sys.stderr.write("Generating consensus sequence for {} individuals\n".
        format(len(sample_list)))

    scaffold_list = []
    if arguments.scaffold is not None:
        scaffold_list = arguments.scaffold.split(',')
    elif arguments.scaffoldfile is not None:
        with open(arguments.scaffold-file,'r') as scaff_file:
            for line in scaff_file:
                scaffold_list.append(line)
    else:
        scaffold_list = None

    #Parse the vcf and write the file
    file_format = arguments.format
    if file_format == 'phylip':
        file_name = arguments.output + '.phylip'
    else:
        file_name = arguments.output + '.fasta'

    fout = open(file_name, 'w')
    if file_format == 'phylip':
        sys.stderr.write("Outputting phylip file\n")
        line_counter = 0
        for sample in sample_list:
            sys.stderr.write("Parsing {}...\n".format(sample))
            variant_alleles = {}
            variant_sites = {}
            variant_alleles, variant_sites, variant_number = parse_vcf(arguments.vcf,sample,
            scaffold_list)

            sys.stderr.write("Writing {}...\n".format(sample))
            if arguments.variants == False:
                if line_counter == 0:
                    fout.write("{}\t{}\n".format(len(sample_list),whole_genome_length))
                    line_counter += 1
                consensus_sequence = ''
                for scaff in reference:
                    if scaffold_list is None or scaff in scaffold_list:
                        new_seq = [x for x in reference[scaff]]
                        for site in range(0,len(variant_sites[scaff])):
                            new_seq[int(variant_sites[scaff][site])-1] = \
                                    variant_alleles[scaff][site]
                        new_seq = ''.join(new_seq)
                        consensus_sequence += new_seq
            elif arguments.variants == True:
                if line_counter == 0:
                    fout.write("{}\t{}\n".format(len(sample_list),variant_number))
                    line_counter += 1
                consensus_sequence = ''
                for scaff in reference:
                    if scaffold_list is None or scaff in scaffold_list:
                        if scaff in variant_alleles:
                            new_seq = ''.join(variant_alleles[scaff])
                            consensus_sequence += new_seq
            fout.write("{}\t{}\n".format(sample,consensus_sequence))
    else:
        for sample in sample_list:
            sys.stderr.write("Parsing {}...\n".format(sample))
            variant_alleles = {}
            variant_sites = {}
            variant_alleles, variant_sites = parse_vcf(arguments.vcf,sample,
            scaffold_list)

            sys.stderr.write("Writing {}...\n".format(sample))

            for scaff in reference:
                if scaffold_list is None or scaff in scaffold_list:
                    new_seq = [x for x in reference[scaff]]
                    for site in range(0,len(variant_sites[scaff])):
                        new_seq[int(variant_sites[scaff][site])-1] = \
                                variant_alleles[scaff][site]
                    new_seq = ''.join(new_seq)
                    fout.write(">{}\n{}\n".format(scaff,new_seq))
    fout.close()

################################################################################


if __name__ == '__main__':
    main()

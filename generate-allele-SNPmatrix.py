#!/usr/bin/env python
'''
Generate a SNP table summarizing allele presence per individual.
The formate here tallies the number of derived alleles present (0,1, or 2).
Missing data is coded as a 9.
'''

import sys
import argparse

def get_args():
    parser = argparse.ArgumentParser(description='Generate SNP frequencies.')
    parser.add_argument(
        '--input',
        required = True,
        help = 'unzipped vcf file'
    )
    parser.add_argument(
        '--output',
        required = False,
        default = 'snp_frequencies.txt',
        help = 'Name for the output file'
    )
    parser.add_argument(
        '--locinames',
        required = False,
        action = 'store_true',
        help = 'generates a list of loci names using the scaffold and pos. \
        No argument required'
    )
    return parser.parse_args()

################################################################################
# Functions

################################################################################

def main():
    arguments = get_args()

    sys.stderr.write("\nParsing VCF file...\n")
    genotype_to_count = { 0:['0/0','./0','0/.'],
                          1:['0/1','1/0','1/.','./1'],
                          2:['1/1'],
                          9:['./.']
    }
    fout = open(arguments.output,'w')
    with open(arguments.input, 'r') as vcf:
        info = []
        header = []
        allele_counts = {}
        num_samples = 0
        num_variants = 0
        for line in vcf:
            if line[0:2] == '##':
                info.append(line)
            elif line.split()[0] == "#CHROM":
                header = line.split()
                num_samples = len(header[9:])
                for id in range(num_samples):
                    id_index = id + 9
                    allele_counts[header[id_index]] = []
            else:
                num_variants+=1
                line = line.split()
                if (num_variants) % 100000 == 0:
                    sys.stderr.write("\nWriting site {}\n".format(num_variants))

                for sample in range(num_samples):
                    sample_index = sample + 9
                    sample_name = header[sample_index]
                    format_field = line[sample_index]
                    genotype = format_field.split(':')[0]
                    for count in genotype_to_count:
                        if genotype in genotype_to_count[count]:
                            allele_counts[sample_name].append(str(count))
    for individual in allele_counts:
        if len(allele_counts[individual]) != num_variants:
            sys.stderr.write('One or more genotypes for ' +
                    '{} are not recognized\n'.format(sample_name))
        else:
            all_counts = '\t'.join(allele_counts[individual])
            fout.write('{}\t{}\n'.format(individual,all_counts))
    fout.close()

    if arguments.locinames == True:
        sys.stderr.write("Writing names of loci to loci_names.txt\n")
        fout_names = open("loci_names.txt", 'w')
        with open(arguments.input, 'r') as vcf_2:
            for line in vcf_2:
                if line[0] != "#":
                    line = line.split()
                    name = []
                    name.append(line[0])
                    name.append(line[1])
                    name_out = ':'.join(name)
                    fout_names.write("{}\n".format(name_out))
        fout_names.close()


################################################################################

if __name__ == '__main__':
    main()

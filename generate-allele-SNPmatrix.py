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

def parse_vcf(input):
    info = []
    header = []
    genos = []
    num_samples = 0
    with open(input,'r') as vcf_file:
        for line in vcf_file:
            if line[0:2] == '##':
                info.append(line)
            elif line.split()[0] == "#CHROM":
                header = line.split()
                num_samples = len(header[9:])
            else:
                line = line.split()
                genos.append(line)
    return info, header, genos, num_samples





################################################################################

def main():
    arguments = get_args()

    sys.stderr.write("\nParsing VCF file...\n")
    info_m, header_m, genos_m, num_samples_m = parse_vcf(arguments.input)

    sys.stderr.write("\nCalculating SNP frequencies...\n")
    genotype_to_count = { 0:['0/0','./0','0/.'],
                          1:['0/1','1/0','1/.','./1'],
                          2:['1/1'],
                          9:['./.']
    }
    fout = open(arguments.output,'w')
    for sample in range(num_samples_m):
        sample_index = sample + 9
        sample_name = header_m[sample_index]
        sample_snp_counts = []
        for snp_line in genos_m:
            format_field = snp_line[sample_index]
            genotype = format_field.split(':')[0]
            for count in genotype_to_count:
                if genotype in genotype_to_count[count]:
                    sample_snp_counts.append(str(count))
        if len(sample_snp_counts) != len(genos_m):
            sys.stderr.write('One or more genotypes for {}'.format(sample_name)+
            'recognized\n')
        else:
            all_counts = '\t'.join(sample_snp_counts)
            fout.write('{}\t{}\n'.format(sample_name,all_counts))
    fout.close()

    if arguments.locinames == True:
        sys.stderr.write("Writing names of loci to loci_names.txt\n")
        fout_names = open("loci_names.txt", 'w')
        for line in genos_m:
            name = []
            name.append(line[0])
            name.append(line[1])
            name_out = ':'.join(name)
            fout_names.write("{}\n".format(name_out))
        fout_names.close()


################################################################################

if __name__ == '__main__':
    main()

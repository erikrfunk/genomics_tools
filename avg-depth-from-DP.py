#!/usr/bin/env python3

'''
Scan a vcf file and calculate the average depth from the DP field for each
variant.

Includes option to filter based on average DP across individuals
'''

import sys
import argparse

def get_args():
    parser = argparse.ArgumentParser(description='Calculate depth per variant')
    parser.add_argument(
        '--vcf',
        required = True,
        help = 'Input vcf'
    )
    parser.add_argument(
        '--output',
        required = False,
        default = 'output.depth',
        help = 'Name of output file'
    )
    parser.add_argument(
        '--field',
        required = False,
        default = 'INFO',
        help = 'Set which field the depth will be read from (FORMAT or INFO) \
        FORMAT option currently not coded!'
    )
    parser.add_argument(
        '--summary',
        required = False,
        action = 'store_true',
        help = 'Set if only a depth summary is desired, and no filtering. \
        If set, no filtering will be done, even if the arguments are set.'
    )
    parser.add_argument(
        '--mindp',
        type = float,
        required = False,
        default = 0,
        help = 'Minimum average depth for a variant'
    )
    parser.add_argument(
        '--maxdp',
        type = float,
        required = False,
        default = -1,
        help = 'Maximum average depth for a variant'
    )
    return parser.parse_args()

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

def calculate_mean_depths(header, genotypes, field, num_samples):
    #Locate the depth value
    field_index = header.index(field)

    mean_depths = []
    if field == 'INFO':
        for variant in genotypes:
            info = variant[field_index].split(";")
            for value in info:
                if value[0:3] == "DP=":
                    depth = value.split("=")[1]
                    mean_depths.append(int(depth) / int(num_samples))
    return mean_depths

###############################################################################

def main():
    arguments = get_args()

    info, header, genos, num_samples = parse_vcf(arguments.vcf)

    sys.stderr.write("\n\nCalculating average depth per variant...\n")

    variant_depth_means = calculate_mean_depths(header, genos, arguments.field,
        num_samples)

    if arguments.summary == True:
        sys.stderr.write("Generating a summary of average depth per variant.\n")
        fout = open(arguments.output, 'w')
        for i in range(len(variant_depth_means)):
            scaffold = genos[i][0]
            pos = genos[i][1]
            fout.write("{}\t{}\t{}\n".format(scaffold, pos,
                variant_depth_means[i]))
        fout.close()
    else:
        if arguments.maxdp == -1:
            maxdp = max(variant_depth_means)
        else:
            maxdp = arguments.maxdp
        sys.stderr.write("Filtering variants shallower than {} and deeper" \
        "than {}.\n".format(arguments.mindp,maxdp))
        fout = open(arguments.output, 'w')
        for line in info:
            fout.write("{}\n".format(line))
        fout.write("{}\n".format("\t".join(header)))

        for i in range(len(variant_depth_means)):
            if (variant_depth_means[i] >= arguments.mindp and
            variant_depth_means[i] <= maxdp):
                genotype = "\t".join(genos[i])
                fout.write("{}\n".format(genotype))
        fout.close()
###############################################################################

if __name__ == '__main__':
    main()

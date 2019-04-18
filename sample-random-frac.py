#!/usr/bin/env python3

'''
Samples a random fraction of variants from a VCF file'''

import sys
import argparse
import random

#Parse the arguments
def get_args():
    parser = argparse.ArgumentParser(description =
    'From a VCF file, sample a random fraction of variants')
    parser.add_argument(
        '--input',
        required = True,
        help = 'Path and name of VCF file'
    )
    parser.add_argument(
        '--fraction',
        required = True,
        help = 'Fraction of total variants to include'
    )
    parser.add_argument(
        '--output',
        required = False,
        default = 'random-variants.vcf'
    )
    return parser.parse_args()


###############################################################################

def main():
    arguments = get_args()

    header = []
    variants = []
    with open(arguments.input,'r') as vcf_file:
        for line in vcf_file:
            if line[0] == "#":
                header.append(line.strip())
            else:
                variants.append(line)

    variant_num = round(len(variants) * float(arguments.fraction))
    sys.stderr.write("Sampling {} variants\n".format(variant_num))

    sampled_variants = [
    variants[i] for i in sorted(random.sample(range(len(variants)),
    variant_num))
    ]

    fout = open(arguments.output, 'w')
    for line in header:
        fout.write("{}\n".format(line))
    for variant in sampled_variants:
        fout.write("{}\n".format(variant))
    fout.close()


###############################################################################

if __name__ == '__main__':
    main()

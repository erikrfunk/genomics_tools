#!/usr/bin/env python
'''
Remove invariant sites from a phylip file.
'''

import sys
import argparse
import numpy

def get_args():
    parser = argparse.ArgumentParser(description='Remove invariant sites from a phylip file.')
    parser.add_argument(
        '--input',
        required = True,
        help = 'phylip file'
    )
    parser.add_argument(
        '--output',
        required = False,
        default = 'trimmed_phylip.phy',
        help = 'Name for the output file'
    )
    return parser.parse_args()

################################################################################

def main():
    arguments = get_args()

    line_num = 0
    seq_length = 0
    sample_num = 0
    seqs = {}
    with open(arguments.input, 'r') as phylip:
        for line in phylip:
            line = line.split()
            if line_num == 0:
                sample_num = line[0]
                seq_length = int(line[1])
                line_num += 1
            else:
                seqs[line[0]] = line[1]

    #Find the invariant sites
    invariant_sites = []
    for pos in range(seq_length):
        if (pos+1) % 100000 == 0:
            sys.stderr.write("\nBase {} of {}\n".format((pos+1),seq_length))
        pos_bases = []
        for sample in seqs:
            pos_bases.append(seqs[sample][pos])
        pos_bases = [x for x in pos_bases if x != 'N']
        if len(set(pos_bases)) == 1:
            invariant_sites.append(pos)

    sys.stderr.write("Removing {} invariant sites.\n".format(len(\
            invariant_sites)))

    fout = open(arguments.output, 'w')
    fout.write("{}\t{}\n".format(sample_num,(seq_length-len(invariant_sites))))
    for sample in seqs:
        seq_list = list(seqs[sample])
        for site in sorted(invariant_sites, reverse = True):
            del seq_list[site]
        seqs[sample] = "".join(seq_list)
        fout.write("{}\t{}\n".format(sample,seqs[sample]))
    fout.close()

################################################################################

if __name__ == '__main__':
    main()

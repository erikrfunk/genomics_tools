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


def del_list_indexes(seq, idxs_to_del):
    trimmed = [i for j, i in enumerate(seq) if j not in idxs_to_del]
    return trimmed

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
            elif line_num == 1:
                if len(line[1]) != seq_length:
                    sys.stderr.write("\nThe sequence length listed in phylip " +
                            "header does not match the length of the of the " +
                            "first sequence.\nYou may want to check this if " +
                            "this file has been used in other analyses.\n\n")
                    seq_length = len(line[1])
                    line_num += 1
                    seqs[line[0]] = line[1]
                else:
                    seqs[line[0]] = line[1]
            else:
                seqs[line[0]] = line[1]

    #Find the invariant sites
    invariant_sites = []
    for pos in range(seq_length):
        if (pos+1) % 100000 == 0:
            sys.stderr.write("\nEvaluating base {} of {}\n".format((pos+1), \
                    seq_length))
        pos_bases = []
        for sample in seqs:
            pos_bases.append(seqs[sample][pos])
        pos_bases = [x for x in pos_bases if x != 'N']
        if len(set(pos_bases)) == 1 or len(set(pos_bases)) == 0:
            invariant_sites.append(pos)

    sys.stderr.write("Removing {} invariant sites.\n".format(len(\
            invariant_sites)))

    fout = open(arguments.output, 'w')
    fout.write("{}\t{}\n".format(sample_num,(seq_length-len(invariant_sites))))

    invariant_sites = set(invariant_sites)
    for sample in seqs:
        sys.stderr.write("Writing {}\n".format(sample))
        seq_list = list(seqs[sample])
        trimmed_seq = del_invars(seq_list,invariant_sites)
        seqs[sample] = "".join(seq_list)
        fout.write("{}\t{}\n".format(sample,seqs[sample]))
    fout.close()

################################################################################

if __name__ == '__main__':
    main()

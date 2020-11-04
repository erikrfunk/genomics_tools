#!/usr/bin/env python3
'''
Take a plink.ld file, and sort values based on their physical bp distance.
Output is a table with column of distances and column of average r2.
'''

import pandas as pd
import numpy as np


def ld_decay(lds,output_file):
    #input_file = "test.ld"
    #output_file = open("chr10_ld_dists.txt",'w')

    snp_dist = [(y-x,z) for x,y,z in zip(lds['BP_A'],lds['BP_B'],lds['R2'])]

    raw_r2_tallies = {}
    for tup in snp_dist:
        if tup[0] not in raw_r2_tallies:
            raw_r2_tallies[tup[0]] = []
            raw_r2_tallies[tup[0]].append(tup[1])
        else:
            raw_r2_tallies[tup[0]].append(tup[1])

    fout = open(output_file, 'w')
    fout.write("Distance\tcount\tr2\n")
    for dist in raw_r2_tallies:
        fout.write("{}\t{}\t{}\n".format(dist, len(raw_r2_tallies[dist]),
        np.mean(raw_r2_tallies[dist])))

    fout.close()

#def ld_by_region(lds,bed,outputfile):
    '''
    bed should be a 3 column file file with "CHR","START","STOP" columns used to
    indicate the ranges the calculation should be broken up into.
    '''

################################################################################
#def main():
     #lds = pd.read_table(input_file,header=0)

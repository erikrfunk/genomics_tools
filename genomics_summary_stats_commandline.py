#!/usr/bin/env python
'''
Calculate some basic genomic diversity and summary statistics.'''

import sys
import argparse

def get_args():
    parser = argparse.ArgumentParser(description='Calculate ')
    parser.add_argument(
        '-i','--input',
        required = True,
        help = 'Fasta file to convert'
    )
    parser.add_argument(
        '-o','--output',
        required = False,
        default='Output.csv'
    )
    parser.add_argument(
        '-p','--popinfo',
        required = False,
        help = 'A file with two columns where the first is sample name and the second is population name'
    )
    parser.add_argument(
        '--pi',
        required = False,
        action = 'store_true',
        help = 'Add flag to calculate sequence diversity for each pop'
    )
    parser.add_argument(
        '--dxy',
        required = False,
        action = 'store_true',
        help = 'Add flag to calculate pairwise sequence divergence across all populations'
    )
    return parser.parse_args()


def FileParse():
    '''Parse the file. Returns a dictionary with sample names as keys and sequences as values'''        
    args = get_args()
    with open (args.input,'r') as f:
        fasta = [line.strip() for line in f.readlines()]   
    sample_names = []
    sequences = []
    for line in fasta:
        if line[0] == ">":
            seq = seq.strip()
            sample_names.append = [seq[1:]]
        else:
            temp_sequence = ''
            while line[0] != ">":
                temp_sequence += line.strip()
            sequences.append(temp_sequence)
    all_seqs = dict(zip(sample_names,sequences))
    return all_seqs

def SetPops():
    '''Use the popinfo file to separate the samples into populations. Returns dictionary.'''
    args = get_args()
    with open(args.popinfo, 'r') as p:
        populations = [line.strip() for line in p.readlines()]
    pops = {}
    for line in populations:
        pops[line.split()[1]] = line.split()[0]
    return pops   

def Pi():
    '''Calculate pairwise nucleotide diversity within a population'''
    all_seqs = FileParse()
    

def Dxy():
    '''Calculate pairwise sequence divergence by averaging pairwise distance across populations'''

def main():
    '''Provide full summary of the file'''
    
    # Start by parsing the file
    # Evaluate the args
    # If pi flagged, calculate pi
    # If dxy flagged, calculate dxy
    
    

if __name__ == '__main__':
    main()

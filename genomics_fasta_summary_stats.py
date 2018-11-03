#!/usr/bin/env python
'''
Calculate some basic genomic diversity and summary statistics.'''

import sys
import argparse
from itertools import combinations

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


def FileParse(Infile):
    '''Parse the file. Returns a dictionary with sample names as keys and sequences as values'''        
    #Add a check to verify all the sequences are the same length
    with open (Infile,'r') as f:
        fasta = [line.strip() for line in f.readlines()]   
    sample_names = []
    sequences = []
    seqcounter = 0
    for line in fasta:
        if line[0] == ">":
            seqcounter += 1
            seq = line.strip()
            sample_names.append(seq[1:])
        else:
            if len(sequences) == (seqcounter - 1):
                sequences.append(line.strip())
            else:
                sequences[(seqcounter - 1)] = sequences[(seqcounter - 1)] + line.strip()
    all_seqs = dict(zip(sample_names,sequences))
    return all_seqs

def SetPops(popsFile):
    '''Use the popinfo file to separate the samples into populations. Returns dictionary.'''
    
    with open(popsFile, 'r') as p:
        populations = [line.strip() for line in p.readlines()]
    pops = {}
    for line in populations:
        if line.split()[1] in pops:
            pops[line.split()[1]].append(line.split()[0])
        else:
            pops[line.split()[1]] = [line.split()[0]]
    return pops   

def Pi(Infile, popsFile, pop):
    '''Calculate pairwise nucleotide diversity within a population'''
    Sequences = FileParse(Infile)
    Populations = SetPops(popsFile)
    print("File contains",len(Sequences),"sequences across",len(Populations),"populations")
    print("and",len(Populations[pop]),"samples in the",pop,"population")
    #print("Sequence are:", Sequences)
    #pops_matrix = [populations[pop],populations[pop]]
    if pop is not None:
        Comparisons = list(combinations(Populations[pop],2))
        #print(Comparisons)
        AlleleFreq = (1/len(Populations[pop]))
        #print(AlleleFreq)
        PiVals = []
        for Pair in Comparisons:
            PairwiseDiffs = []
            #print(Sequences[Pair[0]])
            for Base in range(len(Sequences[Pair[0]])):
                BaseColumn = []
                #print(Sequences[Pair[0]][Base])
                for x in range(len(Pair)):
                    BaseColumn.append(Sequences[Pair[x]][Base])
                PairwiseDiffs.append(len(list(set(BaseColumn))) - 1)
                #print(BaseColumn)
            print(PairwiseDiffs)
            PiVals.append((sum(PairwiseDiffs) / (len(Sequences[Pair[0]]))) * (AlleleFreq**2))
        Pi = (len(Populations[pop]) / (len(Populations[pop])-1)) * (sum(PiVals))
        return Pi
        '''
        for samp1 in pops_matrix[0]:
            for samp2 in pops_matrix[1]:
                print("samples being compared are:", samp1, samp2)
                if samp1 == samp2 or str(samp2 + ' vs ' + samp1) in comps: #avoids duplicate comparisons
                    pass
                else:
                    comps.append(str(samp1 + ' vs ' + samp2))
                    seq_matrix = [sequences[samp1],sequences[samp2]]
                    variability_count = []
                    for basepair in range(len(seq_matrix[1])):
                        base_col = []
                        for seq in seq_matrix:
                            base_col.append(seq[basepair])
                        PairwiseDiffs = (len(set(base_col))-1)/(len(seq_matrix[1]))
                        variability_count.append(PairwiseDiffs*(AlleleFreq**2))
                    pi_vals.append(((sum(variability_count) - len(seq_matrix[0]))/len(seq_matrix[0]))*100)                        
        '''
        average_pi = sum(pi_vals)
    else:
        
        print('Population is variable in ' + str(average_pi) + ' percent of sites')
                
def Dxy(Infile, popsFile, pop1, pop2):
    '''Calculate pairwise sequence divergence by averaging pairwise distance across populations'''


def main():
    '''Provide full summary of the file'''
    
    # Start by parsing the file
    # Evaluate the args
    # If pi flagged, calculate pi
    # If dxy flagged, calculate dxy
    
    

if __name__ == '__main__':
    main()

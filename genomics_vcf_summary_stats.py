#!/usr/bin/env python
'''
Calculate some basic genomic diversity and summary statistics.

Needs to be able to import from VCFparser.py'''

import sys
import argparse
import VCFparser
from VCFparser import Variants
import itertools
from itertools import combinations
'''
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
'''

def Pi(snps_file,pops, pop = None): # Currently broken!! Needs altering to return the right value for pi
    '''Calculate pairwise nucleotide diversity within a population'''
    print("Formatting vcf...")
    Vcf = Variants(snps_file)
    GenomeLength = Vcf.length
    GenotypeLines = Vcf.genotype_convert(genoformat = 'base')
    PopIndex = vcf.index_pops(popsFile = pops)
    print("Calculating average pairwise nucleotide diversity across", Vcf.length, "base pairs")
    Pi = {}
    if pop is not None:
        assert(pop in PopIndex), "Population requested is not listed in the population information file"
        assert(len(PopIndex[pop]) > 1), "Only one individual is provided in the desired population"
        TempPi = []
        Comparisons = list(combinations(PopIndex[pop],2))
        for pair in Comparisons:
            for line in GenotypeLines:
                Genotypes = [line[x] for x in pair]
                Alleles = list(itertools.chain.from_iterable([x.split(sep = '/') for x in Genotypes]))
                Pipq = ((1/(2*len(PopIndex[pop])))**2)()
                TempPi.append(len(set([x for x in Alleles if 'N' not in x])))
        Pi[pop] = ((sum(TempPi)+ GenomeLength)/GenomeLength)
        print('Nucleotide diversity for', pop, 'is', Pi[pop])
        return(Pi[pop])
    else:
        for population in PopIndex:
            if (len(PopIndex[population]) <= 1):
                print('Skipping', population, '- only one individual provided')
            else:
                Pi[population] = []
                Comparisons = list(combinations(PopIndex[population],2))
                for pair in Comparisons:
                    for line in Genotype_lines:
                        Genotypes = [line[x] for x in pair]
                        Alleles = list(itertools.chain.from_iterable([x.split(sep = '/') for x in Genotypes]))
                        Pi[population].append(len(set([x for x in Alleles if 'N' not in x]))-1)
                Pi[population] = ((sum(Pi[population])+ GenomeLength)/GenomeLength)
                print('Nucleotide diversity for', population, 'is', pi[population])
    return(pi)


def Dxy(Infile, popsFile, pop1 = None, pop2 = None):
    '''Calculate pairwise sequence divergence by averaging pairwise distance across populations'''

def Wattersons_theta(snps_file,pops,pop = None):
    '''Watterson's Theta is the number of segregating sites divided by the harmonic number for (n-1) alleles in a given population'''
    #Access the formatted variants using the variants class
    #for every set of length > 1, add one to a segregating sites counter
    #Maybe include an anonymous function to calculate harmonic number
    print("Formatting vcf...")
    Vcf = Variants(snps_file)
    GenomeLength = Vcf.length
    GenotypeLines = Vcf.genotype_convert(genoformat = 'base')
    PopIndex = Vcf.index_pops(popsFile = pops)
    Theta = {}
    if pop is not None:
        SegSites = 0
        assert(pop in PopIndex), "Population requested is not present in the population information file"
        assert(len(PopIndex[pop]) > 1), "Multiple individuals are not present for this population"
        print("Calculating Watterson's Theta for", pop)
        for line in GenotypeLines:
            Genotypes = [line[x] for x in PopIndex[pop]]
            Alleles = list(itertools.chain.from_iterable([x.split(sep = '/') for x in Genotypes]))
            if len(set([x for x in Alleles if 'N' not in x])) > 1: SegSites += 1
        Theta[pop] = (SegSites/GenomeLength)/(sum(map(lambda x: (1/x), list(range(1,len(PopIndex[pop])*2))))) # Multiply by 2 becuase each sample is diploid and a 1 doesn't get subtracted because indexing starts at 0
        print("Watterson's Theta for", pop,":",Theta[pop])
    else:
        print("Calculating Watterson's Theta for", len(PopIndex), "populations")
        for population in PopIndex:
            if (len(PopIndex[population]) <= 1):
                print('Skipping', population, '- only one individual provided')
            else:
                SegSites = 0
                for line in GenotypeLines:
                    Genotypes = [line[x] for x in PopIndex[population]]
                    Alleles = list(itertools.chain.from_iterable([x.split(sep = '/') for x in Genotypes]))
                    if len(set([x for x in Alleles if 'N' not in x])) > 1: SegSites += 1
                Theta[population] = (SegSites/GenomeLength)/(sum(map(lambda x: (1/x), list(range(1,len(PopIndex[population]*2))))))
                print("Watterson's Theta for", population,":",Theta[population])
    return(Theta)





def main():
    '''Provide full summary of the file'''

    # Start by parsing the file
    # Evaluate the args
    # If pi flagged, calculate pi
    # If dxy flagged, calculate dxy



if __name__ == '__main__':
    main()

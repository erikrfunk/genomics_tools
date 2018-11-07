#!usr/bin/env python

import argparse
from VCFparser import Variants

def get_args():
    parser = argparse.ArgumentParser(description='Convert vcf to fasta ')
    parser.add_argument(
        '--input',
        required = True,
        help = 'A vcf file to convert'
    )
    parser.add_argument(
        '--output',
        required = False,
        default='output'
    )
    parser.add_argument(
        '--popinfo',
        required = False,
        default=None,
        help = 'A file with two columns where the first is sample name and the second is population name'
    )
    parser.add_argument(
        '--whichsnp',
        required = False,
        default = 'random',
        help = 'Set how the file will handle diploids: set "pair" to include both the ref and alt alleles in the fasta. Defaul chooses random.'
    )
    parser.add_argument(
        '--individuals',
        required = False,
        default = None,
        help = 'A file with sample names of individuals to include in the fasta')
    return parser.parse_args()

def vcf2fasta(infile, outfile = 'output', popsFile=None,
    whichsnp = 'random', individuals=None, pop=None):
        #Includes two options relating to heterozygotes.
        #'Random' chooses either the ancestral or derived snp randomly
        #'pair' includes both snps, one right after the other with the
        #ancestral listed first
        assert(whichsnp in ['random','pair']), "'whichsnp' needs to be either random or pair. Default is random"
        #Make a dictionary of indices for individuals to be used
        vcf = Variants(infile)
        Samps_to_include = []
        if pop is not None:
            PopIndex = vcf.index_pops(popsFile)
            assert(pop in PopIndex), "Population specified is not in the populations file"
            Samps_to_include = PopIndex[pop]
        #elif individuals is not None:
            #if endswith .txt
            #Open the file with individuals in it and make dict
            #with them and index
            #else set the individual in the same dict
        else:
            print("Writing all samples to fasta...")
            Samps_to_include = vcf.individual_index
        #Strip down the genotypes using the 'genotype_convert' func above
        #Then evaluate which format to create using 'whichsnp'
        #and construct the fasta file accordingly
        formatted_genotypes = vcf.genotype_convert(genoformat = 'base')
        print("There are", str(len(formatted_genotypes)), "variants", sep = ' ')
        reference_index = vcf.header.index('REF')
        if whichsnp == 'pair':
            filename = (outfile + '.fasta')
            with open(filename, 'w+') as fout:
                for sample in Samps_to_include:
                    sequence_header = str('>{0}\n'.format(sample))
                    fout.write(sequence_header)
                    sequence = ''
                    for genotype in formatted_genotypes:
                        temp_geno = genotype[Samps_to_include[sample]].split('/')
                        #If the length of the set is one (homozygous) add the bases
                        if len(set(temp_geno)) == 1:
                            sequence += ''.join(temp_geno)
                        #The the length of the set is two, added the ref snp first, then the alt
                        elif genotype[reference_index] in set(temp_geno):
                            sequence += genotype[reference_index]
                            temp_geno = list(set(temp_geno))
                            temp_geno.remove(genotype[reference_index])
                            sequence += ''.join(temp_geno)
                        else:
                            if len(set(temp_geno)) > 2:
                                print("Genotype > 2. Are you sure all the indels are removed?")
                    fout.write(sequence)
                    fout.write('\n')

        #else: #This should evaluate whichsnp == 'random'

def main():
        arguments = get_args()
        print(arguments)
        vcf2fasta(arguments.input, arguments.output,
                  whichsnp = arguments.whichsnp
        )

if __name__ == '__main__':
    main()

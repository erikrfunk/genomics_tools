#!usr/bin/env python

import re
import gzip

class Variants():
    '''
    Create a vcf class with:

    A modified vcf output that just includes genotypes, 
    A dictionary of individual's index
    A list of indices for individuals in a set of populations
    '''

    def __init__(self, infile = None, genoformat = None, popsfile = None, pop = None, individuals = None):

        self.infile = infile
        self.popsfile = popsfile
        
        with open (infile, 'r') as f:
            vcf = [line.strip() for line in f.readlines()]
        genos = []
        header = []
        length = 0
        for i in vcf:
            text = str(i)
            if re.match('^#[^#]',text):
                header = i.split()
            elif not i.startswith('#'):
                genos.append(i)
            elif re.match('^##contig',text):
                p = re.compile('length=[0-9]+')
                m = p.search(i)
                length += int(m.group()[7:])
                
        self.samplenames = header[9:]
        self.samplenumber = len(self.samplenames)
        self.individual_index = dict(zip(self.samplenames,map(lambda x: header.index(x),self.samplenames)))
        self.rawgenotypes = genos
        self.header = header
        self.length = length
        
        
    def genotype_convert(self, genoformat = None):
        if genoformat not in [None,'raw','numeric','base']:
            print("genotype formate should be genoformat= 'raw', 'numeric' (0/1), or 'base' (A/G)")
        elif genoformat is None:
            return(self.rawgenotypes)
        elif genoformat == 'raw':
            return(self.rawgenotypes)
        elif genoformat == 'numeric':
            new_genos = []
            for i in self.rawgenotypes:
                temp_genos = i.split()
                for field in temp_genos:
                    if temp_genos.index(field) >= (min(list(self.header.index(x) for x in self.samplenames))) - 1:
                        temp_genos[temp_genos.index(field)] = temp_genos[temp_genos.index(field)].split(sep = ':')[0]
                new_genos.append(temp_genos)
            return(new_genos)
        else:
            reference_index = self.header.index('REF')
            alternate_index = self.header.index('ALT')
            new_genos = []
            for i in self.rawgenotypes:
                temp_genos = i.split()
                for field in temp_genos:
                    position = temp_genos.index(field)
                    if position >= (min([self.header.index(x) for x in self.samplenames])):
                        temp_genos[position] = temp_genos[position].split(sep = ':')[0]
                        temp_genos[position] = temp_genos[position].replace('0',temp_genos[reference_index])
                        temp_genos[position] = temp_genos[position].replace('1',temp_genos[alternate_index])
                        temp_genos[position] = temp_genos[position].replace('.','N')
                new_genos.append(temp_genos)
            return(new_genos)

    def index_pops(self, popsFile):
        '''Use the popinfo file or individuals listed in arguments to separate the samples into populations. Returns dictionary.'''
        #Need to add individuals option
        with open(popsFile, 'r') as p:
            populations = [line.strip() for line in p.readlines()]
        pops = {}
        for line in populations:
            if line.split()[1] in pops:
                pops[line.split()[1]].append(self.header.index(line.split()[0]))
            else:
                pops[line.split()[1]] = [self.header.index(line.split()[0])]
        return pops

    def fasta_format(self, outfile = 'output', popsFile=None,
    whichsnp = 'random', individuals=None, pop=None):
        #Includes two options relating to heterozygotes.
        #'Random' chooses either the ancestral or derived snp randomly
        #'pair' includes both snps, one right after the other with the
        #ancestral listed first
        assert(whichsnp in ['random','pair']), "'whichsnp' needs to be either random or pair. Default is random"
        #Make a dictionary of indices for individuals to be used
        Samps_to_include = []
        if pop is not None:
            PopIndex = self.index_pops(popsFile)
            assert(pop in PopIndex), "Population specified is not in the populations file"
            Samps_to_include = PopIndex[pop]
        #elif individuals is not None:
            #if endswith .txt
            #Open the file with individuals in it and make dict
            #with them and index
            #else set the individual in the same dict
        else:
            Samps_to_include = self.individual_index
        print(Samps_to_include)
        #Strip down the genotypes using the 'genotype_convert' func above
        #Then evaluate which format to create using 'whichsnp'
        #and construct the fasta file accordingly
        formatted_genotypes = self.genotype_convert(genoformat = 'base')
        reference_index = self.header.index('REF')
        if whichsnp == 'pair':
            filename = (outfile + '.txt')
            fout = open(filename, 'w+')
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
                    else:
                        sequence += genotype[reference_index]
                        temp_geno.remove(genotype[reference_index])
                        sequence += ''.join(temp_geno)
                fout.write(sequence)
                fout.write('\n')
            fout.close()
        #else: #This should evaluate whichsnp == 'random'

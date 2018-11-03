#!/usr/bin/env python
'''
Convert fasta to nexus.'''

import sys
#import argparse
'''
def get_args():
    parser = argparse.ArgumentParser(description='Convert fasta file to nexus')
    parser.add_argument(
        '--input',
        required = True,
        help = "Fasta file to convert"
    )
    parser.add_argument(
        '--output',
        required = False,
        default='Output.nex'
    )
    return parser.parse_args()

#Get Arguments
#args = get_args()
'''
#------------------------------------------------------------

#If running command line, uncomment this, and comment out below 
#def main():
    #args=get_args()

def main(infile,output="Output.nexus"):    
    with open(infile,'r') as f:
        fasta = [line.strip() for line in f.readlines()]
    
#Parse the file
    sample_num = 0
    sample_names = []
    sequences = []
    for seq in fasta:
        if seq[0] == ">":
            sample_num += 1
            seq = seq.strip()
            sample_names.append(seq[1:])
        else:
            sequences.append(seq)

#Write the taxa block
    nexus = open(output, 'w')
    nexus.write('#nexus\nbegin taxa;\n')    
    temp_string = str('\tdimensions ntax='+str(sample_num)+';\n')
    nexus.write(temp_string)
    nexus.write('\ttaxlabels\n')

    for taxon in sample_names:
        temp_taxon = str('\t\t'+taxon+'\n')
        nexus.write(temp_taxon)
    nexus.write('\t;\nend;\n')

#Write the data block       
    nexus.write('begin data;\n')
    nexus.write(str.format('\tdimentions ntax={0} nchar={1}\n'.format(sample_num,len(max(sequences)))))
    nexus.write('\tformat\tdatatype=dna\tmissing=?\tgap=-;\n\tmatrix\n')
    for i in range(len(sample_names)):
        nexus.write(str.format("\t\t{0}\t{1}\n".format(sample_names[i],sequences[i])))
    nexus.write("\t;\nend;")
    nexus.close()


if __name__=='__main__':
    main(sys.argv[1:])

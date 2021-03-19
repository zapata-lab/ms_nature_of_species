#!/usr/bin/env python

# This script is written by members of the Zapata lab
# Use at your own risk

# UCLA 2019

import os
import sys
import itertools
import argparse
import tempfile

def parse_loci_from_vcf(infile_vcf_file):
    '''Reads filtered vcf_file and returns list of loci'''

    vcf_file = open(infile_vcf_file, "r")

    list_of_loci_in_vcf = []

    for line in vcf_file:
        if line.startswith("locus"):
            list_of_loci_in_vcf.append(line.split()[0])

    return list_of_loci_in_vcf

def parse_original_loci(infile_loci_file):
    '''Reads original loci file from iPyrad and returns list of loci'''

    list_of_locus_numbers = []

    in_loci = open(infile_loci_file, "r")

    for line in in_loci:
      if line.startswith("//"):
        line = line.split()
        list_of_locus_numbers.append("locus_"+ line[-1].strip("|*-"))

    return list_of_locus_numbers
    
def parse_names_to_keep(infile_vcf_file):
	'''Reads line from vcf_file that holds names of samples'''
	
	vcf_file = open(infile_vcf_file, "r")
	
	for line in vcf_file:
		if line.startswith("#CHROM"):
			list_of_names_to_keep = line.split()[9:]
	
	return list_of_names_to_keep

def write_new_gphocs(infile_loci_file, infile_vcf_file, list_of_loci_wanted, list_of_loci_numbers_from_original_file, list_of_kept_names):
    '''Uses loci from filtered vcf file and sequences from the original loci file \
    to create a gphocs file and write the new filtered gphocs file to the current directory.'''

    original_loci_file = open(infile_loci_file, "r")
    filename = os.path.basename(infile_vcf_file).split(".")[0]
    outfile =  os.path.join(os.getcwd(), filename + "_handmade.gphocs")

    loci = original_loci_file.read().strip().split("//")[:-1]

    with tempfile.TemporaryFile() as holder:
        for locus_info, specific_locus_number in itertools.izip(loci, list_of_loci_numbers_from_original_file):
            names = [line.split()[0] for line in locus_info.strip().split("\n")]
            sequences = [line.split()[-1] for line in locus_info.strip().split("\n")]
            sequences = [seq.replace("n","").replace('-','N') for seq in sequences]
            for name,sequence in zip(names,sequences):
                if name not in list_of_kept_names:
                    i = names.index(name)
                    names.pop(i)
                    sequences.pop(i)
            #print names
            #print sequences
            #print len(sequences)
	    if len(sequences) > 0:
            	sequence_length = len(sequences[0])
	    	longname = max(map(len,names))+ 4
                print >> holder, "\n", specific_locus_number, len(names), str( sequence_length )
                for name,sequence in zip(names,sequences):
                    if len(sequence) > 20:
                        print >> holder, name+" "*(longname-len(name))+sequence
        holder.seek(0)

        with open(outfile, "w") as f:
            print >> f, len(list_of_loci_wanted)

        for line in holder:
          if line.startswith("locus"):
            line_cont = line.rstrip().split(" ")
            if line_cont[0] not in list_of_loci_wanted:
                line = next(holder)
            else:
                with open(outfile, "a") as f:
                    print >> f, "\n", line.rstrip()
                    line = next(holder)
                    while not line.startswith("\n"):
                        print >> f, line.rstrip()
                        line = next(holder)
	
	with open(outfile, "a") as f:
        	print >> f, "\n"


if __name__ == "__main__":
    parser = argparse.ArgumentParser( description = 'This script filters an original loci file (from iPyrad), based on loci found in a filtered vcf file, and generates a file in the gphocs format.')
    parser.add_argument('--infile_vcf_file', help = """Final vcf file from filtering based on whatever filtering scheme you use in vcftools. Provide full path if running this script from a different directory.""")
    parser.add_argument('--infile_loci_file', help = """Original loci file that has all loci (from an iPyrad processing pipeline). Provide full path if running this script from a different directory.""")
    args = parser.parse_args()

    listofdesiredloci = parse_loci_from_vcf(args.infile_vcf_file)
    listoforiginallocinumbers = parse_original_loci(args.infile_loci_file)
    listofkeptnames = parse_names_to_keep(args.infile_vcf_file)
    write_new_gphocs(args.infile_loci_file, args.infile_vcf_file, listofdesiredloci, listoforiginallocinumbers, listofkeptnames)

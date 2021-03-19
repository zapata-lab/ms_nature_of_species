#!/usr/bin/env python

# This script is written by members of the Zapata lab
# Use at your own risk

# UCLA 2019

import argparse
import os

def write_partition_file(infile_gphocs_file):
	
	gphocs_file = open(infile_gphocs_file, "r")
	
	locus_info = []
	
	for line in gphocs_file:
		if line.startswith("locus"):
			locus_info.append(line.split())
	base_pair_num = 0
	total_length = []
	filename = os.path.basename(infile_gphocs_file).split(".")[0]
	outfile =  os.path.join(os.getcwd(), filename + ".par")
	with open(outfile, "w") as f:
		for i in locus_info:
			print >> f, "DNA, ", i[0], "=", int(base_pair_num) + 1, "-", int(base_pair_num) + int(i[2])
			total_length.append(int(i[2]))
			base_pair_num = sum(total_length)



def write_locus_files(infile_gphocs_file):
	
	gphocs_file = open(infile_gphocs_file, "r")
	iterofile = iter(gphocs_file)

	for line in iterofile:
		if line.startswith("locus"):
			line_cont = line.split(" ")
			filename = line_cont[0]
			outfile = os.path.join(os.getcwd(), filename + ".phy")
			with open(outfile, "w") as f:
				print >> f, line_cont[1], line_cont[2]
				line = next(iterofile)
				line = line.rstrip()
				while not line.startswith('\n'):
					print >> f, line.rstrip()
					line = next(iterofile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser( description = 'This script creates a partition file from the loci listed in a gphocs file in the format')
    parser.add_argument('--infile_gphocs_file', help = """Gphocs file. Provide the full path to the file if running from a different directory.""")
    args = parser.parse_args()

    partitionfile = write_partition_file(args.infile_gphocs_file)
    locusfiles = write_locus_files(args.infile_gphocs_file)
	

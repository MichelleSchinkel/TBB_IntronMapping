#!/usr/bin/env python3

import sys; import getopt
from Bio import SeqIO

def usage():
    print("\n\npath_to_script.py [-h] [options] <inputfile.tsv>\n")
    print("<inputfile.tsv> Needs to be a table with OG's in column 1 and eukarya-4-identifiers in column 2\n")
    print("Options:")
    print("-o: output file (default: output.fasta)")
    print("-e: eukarya-4-database file (default:/home/julian/julian2/snel-clan-genomes/eukarya/eukarya_4.0.1/eukarya_4.0.1_proteomes_lt.fa)")
    
optlist, args = getopt.getopt(sys.argv[1], 'o:e')

opts = {}
for k,v in optlist:
    if k == '-h':
        usage(); sys.exit()
    else:
        opts[k] = v  

OG_doc = open(sys.argv[1], "r")

# Defining the output file 
if '-o' in opts.keys():
    output_path = opts["-o"]
else:
    output_path = "output.fasta"
output_doc = open(output_path , "w")

# Defining where to find the eukarya-4-database
if "-e" in opts.keys():
    eukarya_path = opts["-e"]
else:
    eykarya_path = "/home/julian/julian2/snel-clan-genomes/eukarya/eukarya_4.0.1/eukarya_4.0.1_proteomes_lt.fa"
eukarya_doc = open(eukarya_path , "r")


tubulin_dic = {}

for line in OG_doc:
    line = line.strip()
    words = line.split("\t")
    tubulin_dic[words[1]] = str(words[0])
    

for record in SeqIO.parse(eukarya_doc, 'fasta'):
    for key in tubulin_dic.keys():
        if key in record.id:
            print(">" + tubulin_dic[key] + "_" + key, file=output_doc)
            print(record.seq + "\n", file=output_doc)

OG_doc.close()
eukarya_doc.close()
output_doc.close()



#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 11:57:53 2019

@author: michelle
"""
from Bio import SeqIO
import os

lecas = open("/home/michelle/michelle2/pfam/LECA-OGs/lecas_all_seqs.tsv", "r")



for line in lecas:
    if "Support" in line:
        pass
    else:
        line = line.strip()
        line = line.split("\t")
        pfam_number = line[0]
        OG_donation = line[1]
        stop = OG_donation.rfind(".")
        try:
            os.mkdir("/home/michelle/michelle2/pfam/LECA-OGs/" + pfam_number)
        except OSError:
            pass
        output_doc = open("/home/michelle/michelle2/pfam/LECA-OGs/" + pfam_number + "/" + pfam_number + "_" + OG_donation + ".fa", "a")
        eukarya4 = line[3]
        eukarya4 = eukarya4.split(",")
        for record in SeqIO.parse("/home/michelle/michelle2/pfam/sequences/" + pfam_number + "_2_domains_numbered.fa", "fasta"):
            for identifier in eukarya4:
                if identifier in record.id:
                    print(">" + OG_donation + "_" + identifier, file=output_doc)
                    print(record.seq, file=output_doc)
    
    
os.mkdir
    
    
    
    



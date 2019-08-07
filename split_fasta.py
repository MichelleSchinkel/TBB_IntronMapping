#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 12:54:50 2019

@author: michelle
"""

numbers = list(range(401, 563))

for number in numbers:
    number = str(number)
    fasta_file = open("/home/michelle/Documents/project/data/sverdlov_kogs/all_orth_groups/group" + number + "/group" +number + ".fasta", "r")


    for line in fasta_file:
        line = line.strip()
        if line.startswith(">") == False:
            print(line, file=kog_file)
            kog_file.close()
        if line.startswith(">"):
            end = line.rfind("_")
            kog = line[1:end]
            kog_file = open("/home/michelle/Documents/project/data/sverdlov_kogs/all_orth_groups/group" +number+"/" + kog + ".fasta", "a+")
            print(line, file=kog_file)


               
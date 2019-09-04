#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 16:40:13 2019

@author: michelle
"""

import os

non_trusted = open("/home/michelle/michelle2/pfam/LECA-OGs/non_trusted/total_list.txt", "r")

for line in non_trusted:
    line = line.strip()
    words = line.split("\t")
    pfam = words[0]
    OG = words[1]
    try:
        os.mkdir("/home/michelle/michelle2/pfam/LECA-OGs/non_trusted/" + pfam)
    except OSError:
        pass

    numbers = list(range(1, 99))

    for number in numbers:
        number = str(number)

        try:
            os.rename("/home/michelle/michelle2/pfam/LECA-OGs/non_duplicated/" + pfam + "/" + pfam + "_" + OG + "." + number + ".fa" , "/home/michelle/michelle2/pfam/LECA-OGs/non_trusted/" + pfam + "/" + pfam + "_" + OG + "." + number + ".fa" )  
            os.rename("/home/michelle/michelle2/pfam/LECA-OGs/non_duplicated/" + pfam + "/" + pfam + "_" + OG + "." + number + ".aln" , "/home/michelle/michelle2/pfam/LECA-OGs/non_trusted/" + pfam + "/" + pfam + "_" + OG + "." + number + ".aln" )
        except OSError:
            pass

non_trusted.close()
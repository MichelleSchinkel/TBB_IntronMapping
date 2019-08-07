#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 14:35:03 2019

@author: michelle
"""

#tubulin OGs: tubulinx, eukarya-4-identifier
#KOG OGs: KOG, eukarya-4-identifier
tubulin_OGs = open("/home/michelle/Documents/project/data/tubulin/tubulin_LECA-OGs.tsv", "r")
KOG_OGs = open("/home/michelle/Documents/project/data/sverdlov_kogs/all_orth_groups/group181/group181.txt", "r")
tub_KOG = open("/home/michelle/Documents/project/data/sverdlov_kogs/tub_KOG.txt", "w")

tub_dic = {}

for line in tubulin_OGs:
    line = line.strip()
    words = line.split()
    tub_dic[words[1]] = str(words[0])
#    print(tub_dic)

for line in KOG_OGs:
    line = line.strip()
    words2 = line.split()
    KOG = words2[0]
    eukarya = words2[1]
    for key in tub_dic.keys():
        if key == eukarya:
            print(KOG, "\t", tub_dic[key], "\t", eukarya, file=tub_KOG)

tubulin_OGs.close()
KOG_OGs.close()
tub_KOG.close()



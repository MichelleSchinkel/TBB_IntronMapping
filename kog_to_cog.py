#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 10:37:30 2019

@author: michelle
"""






kog_to_cog = open("/home/michelle/Documents/project/data/sverdlov_kogs/kog_jolien_cut.txt", "r")
kog_to_cog_dic = {}
words = ""

for line in kog_to_cog:
    line = line.strip()
    words = line.split("\t")
    kog = words[0]
    cog = words[1]
    kog_to_cog_dic[kog] = [cog]
    words = ""
#print(kog_to_cog_dic)

kog_to_cog.close()


kog_groups = open("/home/michelle/Documents/project/data/sverdlov_kogs/kog_groups_suppl.txt", "r")
cog_groups = open("/home/michelle/Documents/project/data/sverdlov_kogs/cog_groups.txt", "w+")

for line in kog_groups:
    line = line.strip()
    if line.startswith(">"):
        print(line, file=cog_groups)
    else:
        words = line.split(" ")
        for word in words:
            if word in kog_to_cog_dic.keys():
                if word not in cog_groups:
                    print (kog_to_cog_dic[word], file=cog_groups)

#print(cog_groups)

kog_groups.close()
cog_groups.close()   
    
#cog_groups = open("/home/michelle/Documents/project/data/sverdlov_kogs/cog_groups.txt", "r")
#cog_groups_uniq = open("/home/michelle/Documents/project/data/sverdlov_kogs/cog_groups_uniq.txt", "w")
#             
#for line in cog_groups:
#    word = line.split()
#    if word not in cog_groups:
#        print(word, file=cog_groups_uniq)              
#
#cog_groups.close()
#cog_groups_uniq.close()

            
        

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 12:42:01 2019

@author: Michelle Schinkel
"""

#####################################
# Create a dictionary of the eggnog annotations wiht kog as key and eukarya4 identifiers as values.
eggnog = open("/home/michelle/michelle2/eukarya/eukarya/annotations/eggnog/EggNOGgroups.txt", "r")

eggnog_dic = {}

for line in eggnog:
    line = line.strip()
    words = line.split(" ")
    kog_eunog = words[0]
    kog = kog_eunog[0:7]
    eggnog_dic[kog] = words[1:]

eggnog.close()

#####################################
# Couple the kog's to the values of the eggnog dictionary to retrieve the eukarya4 identifiers per group
kog_groups = open("/home/michelle/Documents/project/data/sverdlov_kogs/kog_groups_suppl.txt", "r")
eukarya4_groups = open("/home/michelle/Documents/project/data/sverdlov_kogs/eukarya4_groups.txt", "w")

for line in kog_groups:
    line = line.strip()
    if line.startswith(">"):
        print(line, file=eukarya4_groups)
    else:
        words = line.split(" ")
        for word in words:
            if word in eggnog_dic.keys():
                for value in eggnog_dic[word]:
                    print (word + "\t" + str(value), file=eukarya4_groups)

kog_groups.close()
eukarya4_groups.close()

######################################
## Trimming the document by removing trailing and unnecessary characters.
#eukarya = open("/home/michelle/Documents/project/data/sverdlov_kogs/eukarya4_groups.txt", "r")
#eukarya4_trimmed = open("/home/michelle/Documents/project/data/sverdlov_kogs/eukarya4_groups_trimmed.txt", "w")
#
#for line in eukarya:
#    line = line.strip()
#    line = line.replace("[","")
#    line = line.replace("]","")
#    line = line.replace(",","")
#    line = line.replace("'","")
#    print(line, file=eukarya4_trimmed)
#
#eukarya4_groups.close()
#eukarya4_trimmed.close()


#######################################

eukarya4 = open("/home/michelle/Documents/project/data/sverdlov_kogs/eukarya4_groups.txt", "r")

for line in eukarya4:
    if line.startswith(">"):
        line = line.strip()
        file = open("/home/michelle/Documents/project/data/sverdlov_kogs/all_orth_groups/"+str(line[1:])+".txt" , "w+")
    else:
        line = line.strip()
        print(line, file=file)

file.close()
eukarya4.close()

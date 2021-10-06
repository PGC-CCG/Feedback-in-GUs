#!/usr/bin/env python
# coding: utf-8

###################################################
## Create Regulons                               ##
## Author: Georgette Femerling                   ##
## Version: v1                                   ##
## Description: Script that creates simple       ##
##          Regulons from a tf-gen network       ##
## Input: [1] Input Network                      ##
##        [2] Output name                        ##
## Output: Regulon Network with the format:      ##
## TF\tGene1//gene2//....//                      ##
###################################################

#import modules
from __future__ import absolute_import, division
import re
import networkx as nx
import sys
import os
from datetime import date

def usage():
    """help message function"""

    print """ 
Create_Regulons.py : Creates general regulons from a tab-separated tf-gen regulatory network file. Output regulon is a list of genes separated by //.

python Create_Regulons.py [-h] <Input Network> <Output name>

    -h          print this message
    <Input Network>  required. tf-gen network to convert. Must be tab separated.
    <Output name> required. Name of output file

        """

if len(sys.argv) < 2:
    usage(); sys.exit("Create_Regulons.py: Required arguments missing\n")

if "-h" in sys.argv:
    usage(); sys.exit()

if len(sys.argv) < 3:
    usage(); sys.exit("Create_Regulons.py: Required arguments missing\n")

#Get arguments
input_ntw = sys.argv[1]
output_ntw = sys.argv[2]

Ecoli = nx.DiGraph()
with open (input_ntw) as archivo: 
    for line in archivo:
        line=line.rstrip()
        if not re.search('#',line):
            TF=line.split()[0] #Se hace un split con tab y se toma el primer elemento - el TF
            Gene=line.split()[1] #tomo el segundo elemento del split - el gen Regulado
            Ecoli.add_edge(TF, Gene)

TF_list=[]
for TF in Ecoli.out_degree():
    if TF[1]>0:
        TF_list.append(TF[0])

#crear los regulones
Regulones={}
for TF in TF_list:
    Regulones[TF] = Ecoli[TF].keys()

#Escribir el archivo
newfile=open(output_ntw+'.txt','w')
newfile.write("# From Create_Regulons.py \n")
newfile.write("# Date: "+str(date.today())+"\n")
newfile.write("# Input network: "+input_ntw+"\n")
for TF in Regulones.keys():
    newfile.write(TF+'\t'+"//".join(Regulones[TF])+"//"+'\n')
    
print(".........Done..........\n")


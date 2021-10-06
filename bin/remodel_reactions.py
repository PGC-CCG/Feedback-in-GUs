#!/usr/bin/env python
# coding: utf-8
# Import Modules
import re
import sys
import os
from datetime import date

###################################################
## Remodel Reactions                             ##
## Author: Georgette Femerling                   ##
## Version: v1                                   ##
## Date: 8-04-2020                               ##
## Description: Script that takes the output     ##
## files from the GU and generates one file with ##
## the reactions in linear format                ##
## Input: [1] GUs_secondary_reactions dir        ##
##        [2] complete path to Output_dir        ##
##            to create                          ##
## Output: One file per GU with lineal reactions ##
###################################################

### Example of run python ./remodel_reactions.py ../GUs_fromgroup/PhysicalNt_groups2/datasets/GUs_secondary_reactions/ ../GUs_fromgroup/PhysicalNt_groups2/datasets/Remodelled_GUs

"""Remodel Reactions: Takes the output files from the GU assembly and generates one file with the GU reactions in linear format """

def usage():
    """help message function"""

    print """ 
remodel_reactions.py : reads GU assembly output files and generates one file per GU with the reactions in linear format.

python remodel_reactions.py [-h] <GU_secondary_reactions dir> <Output dir>

    -h          print this message
    <GU_secondary_reactions dir>  required. GU secondary reactions directory
    <Output dir> required. Name of output directory

        """

if "-h" in sys.argv:
    usage(); sys.exit()

if len(sys.argv) < 3:
    usage(); sys.exit("remodel_reactions.py: Required arguments missing\n")

# Get input directory
input_dir = sys.argv[1] #input_dir = "../GUs_fromgroup/PhysicalNt_groups2/datasets/GUs_secondary_reactions/"
output_dir = sys.argv[2] #output_dir = "../GUs_fromgroup/PhysicalNt_groups2/datasets/Remodelled_GUs"
try:
    os.mkdir(output_dir) #Creates output dir
except:
    print str(output_dir)+" already exists\n"

# Main Code
for GU in os.listdir(input_dir):
    if GU == ".DS_Store":
        continue
    print "........{}........\n".format(GU)
    rxs = {}
    # Read reactant_product file and arrange the values in a dictionary of reactants and products per reaction
    with open (os.path.join(input_dir,GU+"/reactants_products.txt")) as rctprd:
        for line in rctprd:
            rx = re.match("^(\w+)\t(.+)\t(\w+)",line) #parse line
            if rx.group(1) in rxs.keys(): # if re# in the dictionary
                if rx.group(3) == 'reactant': #add reactant
                    rxs[rx.group(1)]['reactants'].append(rx.group(2))
                elif rx.group(3) == 'product': #add product     
                    if 'products' in rxs[rx.group(1)].keys():
                        rxs[rx.group(1)]['products'].append(rx.group(2))
                    else:
                        rxs[rx.group(1)]['products'] = [rx.group(2)]
            else:
                rxs[rx.group(1)] = {} #create it
                if rx.group(3) == 'reactant': #add reactant
                    rxs[rx.group(1)]['reactants'] = [rx.group(2)]

    # Read reactions file to get the type of reaction of each reaction
    with open (os.path.join(input_dir,GU+"/reactions.txt")) as rxsfile:
        for line in rxsfile:
            rx = re.match("^(\w+)\t(\w+)",line) #parse line
            rxs[rx.group(1)]['type'] = rx.group(2) #add type of reaction

    # Read modifications file to get the subtyope of reaction
    with open (os.path.join(input_dir,GU+"/modification.txt")) as mod:
        for line in mod:
            rx = re.match("^(\w+)\t(\w+)\t?([\w\-\_\|\(\)]+)?",line) #parse line
            rxs[rx.group(1)]['subtype'] = rx.group(2) #Add subtype
            if rx.group(3) not in rxs[rx.group(1)]['reactants']: #Add extra reactant if missing
                rxs[rx.group(1)]['reactants'].append(rx.group(3))

    # Read complexes file to get the complex formation reactions           
    with open (os.path.join(input_dir,GU+"/complexes.txt")) as comp:
        for line in comp:
            cs = re.match("^(\w+)\t([\w\-]+)\t(\w+)",line) #parse line
            for rx in rxs:
                if 'products' in rxs[rx].keys():
                    if cs.group(2) in rxs[rx]["products"]: #Verify if the complex is a product
                         rxs[rx]["subtype"] = "COMPLEX" #Add subtype
                else:
                    rxs[rx]["products"] = ""

    # Write modelled reactions file
    newfile = open(os.path.join(output_dir,GU+"_Modelled.txt"),"w")
    # Write header
    newfile.write("# From remodel_reactions.py \n")
    newfile.write("# Date: "+str(date.today())+"\n")
    newfile.write("# Reactants\tProducts\tReaction_Type\n")
    for ren in rxs:
        reactants = rxs[ren]['reactants']
        products = rxs[ren]['products']
        if rxs[ren]['type'] == "TRANSCRIPTION":
            rxtype = "Transcription"
        elif rxs[ren]['type'] == "TRANSLATION":
            rxtype = "Translation"
        elif rxs[ren]['type'] == "TRANSPORT":
            rxtype = "Transport"
        elif rxs[ren]['type'] == "STATE_TRANSITION":
            if 'subtype' in rxs[ren].keys():
                if rxs[ren]['subtype'] == "COMPLEX":
                    rxtype = "Complex"   
                else:
                    rxtype = "Catalysis"
        elif rxs[ren]['type'] == "SUPER_REACTION":
            rxtype = "Catalysis"
        newfile.write(",".join(reactants)+"\t"+",".join(products)+"\t"+rxtype+"\n")
    newfile.close()

print "Done\n"

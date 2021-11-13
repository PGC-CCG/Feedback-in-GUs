#!/usr/bin/env python3

# Import Modules
import re
from collections import defaultdict
from datetime import date
import os
import sys

conformationf = "./gu_library/Conformaciones_TF_efector.tab.txt"
# GU_dir = "/home/gfemer/GUs/GUs_Regulon10.7_40821"
# conformationf = sys.argv[1]
GU_dir = sys.argv[1]
outputdir = sys.argv[2]

# Read TF conformation table
conformation_dic = defaultdict(set)
conformation_dic2 = defaultdict(set)
effector_dic = defaultdict(int)            
with open(conformationf) as conf:
    for line in conf:
        # Skip first line - header
        if(re.match("^#",line) or re.match("^TRAN_FACTOR",line)):
            continue
        line = line.strip("\n").strip("\r").split("\t")
        TF = line[0]
        effector_dic[TF]+=1
        # Remove weird symbols from effector name
        Effector_name = re.sub(r"[\"&;|]|(</?.+?>)","",line[7])
        conformation_dic[TF].add(Effector_name)
        conformation_dic2[TF].add(Effector_name)
        # Get any synonyms, also removing symbols
        if len(line[9]) > 0:
            synonyms = re.sub(r"(\'),(\d\')",r"\1_\2",re.sub(r"[\"&;|]|(</?.+?>)","",line[9])).split(",")
            # This is to account for the comma in conformations with carbon numbers such as 5',3'-cAMP 
            synonyms = [re.sub(r'_', ',', i) for i in synonyms]
            conformation_dic[TF].update(synonyms)

# Actually Calculating the feedback
def feedback(dir,conf_dic):
    TFs_in_net = set()
    feedback_dic = defaultdict(lambda: defaultdict(int)) #Dictionary to keep track of how many loops are found
    matches_dic = defaultdict(lambda: defaultdict(set)) #To keep track of possible matches if not exact
    for TF in os.listdir(os.path.join(dir+"/GUs_secondary_reactions")):
        TFs_in_net.add(TF) #Added TF to the set
        if TF not in conf_dic.keys():
            continue
        with open(os.path.join(dir+"/GUs_secondary_reactions/"+TF+"/reactants_products.txt")) as rp:
            found = [] #Flag and to not repeat same object
            matches = []
            for line in rp:
                transp=0
                obj = line.strip().split("\t")[1]
                if re.search("_Ext",obj):
                    transp=1
                obj = re.sub(r"[\"&;|]|(</?.>)|_Ext",'',obj)
                # Exact match
                if (((obj in conf_dic[TF]) and (obj not in found)) or ((obj in conf_dic[TF]) and (transp == 1) and (feedback_dic[TF][obj] == 0))):
                    found.append(obj)
                    feedback_dic[TF][obj] = transp
        # If no feedback found, search for closest compound           
        if len(found) == 0:
            with open(os.path.join(dir+"/GUs_secondary_reactions/"+TF+"/reactants_products.txt")) as rp:
                # print(TF)
                for line in rp:
                    transp = 0
                    obj = line.strip().split("\t")[1]
                    if re.search("_Ext",obj):
                        transp=1
                    obj = re.sub(r"[\"&;|]|(</?.>)|_Ext",'',obj)
                    # If an effector contains the object as a string
                    obj2 = re.sub(r"L-?|D-?|alpha[,-]|beta[,-]|gamma[,-]|omega[,-]|delta[,-]|","",obj)
                    matches = [cmp for cmp in conf_dic[TF] if obj2 == re.sub(r"L-?|D-?|alpha[,-]|beta[,-]|gamma[,-]|omega[,-]|delta[,-]|","",cmp)]
                    if len(matches) > 0:
                        for effector in matches:
                            matches_dic[TF][effector].add(obj)
                            if ((obj not in feedback_dic[TF].keys()) or (transp == 1 and (obj in feedback_dic[TF].keys()) and feedback_dic[TF][obj] == 0)):
                                feedback_dic[TF][obj] = transp
                    # If an object contains the effector as a string
                    # matches2 = [cmp for cmp in conf_dic[TF] if cmp in obj]
                    # if len(matches2) > 0:
                    #     for effector in matches2:
                    #         matches_dic[TF][effector].add(obj)
    # Return dictionary with TFs with exact feedback and dictionary of possible feedback
    return feedback_dic,TFs_in_net,matches_dic

feed_dic,TF_list,matches = feedback(GU_dir,conformation_dic)

with open(os.path.join(outputdir+"/TFs_with_feedback.txt"),"w") as fd:
    fd.write("# TF\tEfector(s)\tNo. of Efectors\tMatches\tNo. of matches\tTransported Effectors\n")
    for TF in feed_dic:
        fd.write(str(TF)+"\t"+",".join(conformation_dic2[TF])+"\t"+str(len(conformation_dic2[TF]))+"\t"+",".join(feed_dic[TF].keys())+"\t"+str(len(feed_dic[TF].keys()))+"\t"+",".join([key for key in feed_dic[TF].keys() if feed_dic[TF][key]==1])+"\n")

with open(os.path.join(outputdir+"/TFs_no_feedback.txt"),"w") as nf:
    nf.write("# TF\tEfector\n")
    TFs=set(conformation_dic.keys()).intersection(TF_list)
    # nf_set=set(conformation_dic.keys()).difference(set(feed_dic.keys()))
    nf_set=TFs.difference(set(feed_dic.keys()))
    for TF in nf_set:
        nf.write(str(TF)+"\t"+list(conformation_dic[TF])[0]+"\n")

with open(os.path.join(outputdir+"/TFs_feedback_possible_matches.txt"),"w") as pm:
    pm.write("# TF\tEf.synonym\tmatch\n")
    for TF in matches:
        for syn in matches[TF]:
            pm.write(str(TF)+"\t"+str(syn)+"\t"+",".join(matches[TF][syn])+"\n")

print("Number of TFs with GU: "+str(len(TF_list)))
print("Number of TFs with known effector: "+str(len(effector_dic.keys())))
print("Number of TFs with feedback: "+str(len(feed_dic.keys())))
print("Percentage of feedback in the network is: "+str(len(feed_dic.keys()) / len(effector_dic.keys())))
print("\nDone")
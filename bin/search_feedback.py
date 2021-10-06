
#!/usr/bin/env python3

# Import Modules
import re
from collections import defaultdict
from datetime import date
import os

conformationf = "/home/gfemer/GUs/gu_library/Conformaciones_TF_efector.tab.txt"
GU_dir = "/home/gfemer/GUs/GUs_Regulon10.7_40821"

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
    feedback_dic = defaultdict(list) #Dictionary to keep track of how many loops are found
    matches_dic = defaultdict(lambda: defaultdict(set)) #To keep track of possible matches if not exact
    for TF in os.listdir(os.path.join(dir+"/GUs_secondary_reactions")):
        if TF not in conf_dic.keys():
            continue
        with open(os.path.join(dir+"/GUs_secondary_reactions/"+TF+"/reactants_products.txt")) as rp:
            found = [] #Flag and to not repeat same object
            matches = []
            for line in rp:
                obj = line.strip().split("\t")[1]
                obj = re.sub(r"[\"&;|]|(</?.>)|_Ext",'',obj)
                # Exact match
                if (obj in conf_dic[TF]) and (obj not in found):
                    found.append(obj)
                    feedback_dic[TF].append(obj)
        # If no feedback found, search for closest compound           
        if len(found) == 0:
            with open(os.path.join(dir+"/GUs_secondary_reactions/"+TF+"/reactants_products.txt")) as rp:
                print(TF)
                for line in rp:
                    obj = line.strip().split("\t")[1]
                    obj = re.sub(r"[\"&;|]|(</?.>)|_Ext",'',obj)
                    # If an effector contains the object as a string
                    # matches = [cmp for cmp in conf_dic[TF] if obj in cmp]
                    obj2 = re.sub(r"L-?|D-?|alpha[,-]|beta[,-]|gamma[,-]|omega[,-]|delta[,-]|","",obj)
                    matches = [cmp for cmp in conf_dic[TF] if obj2 == re.sub(r"L-?|D-?|alpha[,-]|beta[,-]|gamma[,-]|omega[,-]|delta[,-]|","",cmp)]
                    if len(matches) > 0:
                        for effector in matches:
                            matches_dic[TF][effector].add(obj)
                            if obj not in feedback_dic[TF]:
                                print(effector,obj)
                                feedback_dic[TF].append(obj)
                    # If an object contains the effector as a string
                    # matches2 = [cmp for cmp in conf_dic[TF] if cmp in obj]
                    # if len(matches2) > 0:
                    #     for effector in matches2:
                    #         matches_dic[TF][effector].add(obj)
    # Return dictionary with TFs with exact feedback and dictionary of possible feedback
    return feedback_dic,matches_dic

feed_dic,matches = feedback(GU_dir,conformation_dic)

with open("TFs_with_feedback.txt","w") as fd:
    fd.write("# TF\tEfector(s)\tNo. of Efectors\tMatches\tNo. of matches\n")
    for TF in feed_dic:
        fd.write(str(TF)+"\t"+",".join(conformation_dic2[TF])+"\t"+str(len(conformation_dic2[TF]))+"\t"+",".join(feed_dic[TF])+"\t"+str(len(feed_dic[TF]))+"\n")

with open("TFs_no_feedback.txt","w") as nf:
    nf_set=set(conformation_dic.keys()).difference(set(feed_dic.keys()))
    for TF in nf_set:
        nf.write(str(TF)+"\t"+list(conformation_dic[TF])[0]+"\n")

with open("TFs_no_feedback_possible_matches.txt","w") as pm:
    pm.write("# TF\tEf.synonym\tmatch\n")
    for TF in matches:
        for syn in matches[TF]:
            pm.write(str(TF)+"\t"+str(syn)+"\t"+",".join(matches[TF][syn])+"\n")
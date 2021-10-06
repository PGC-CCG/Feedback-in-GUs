#!/usr/bin/env python3

# Import Modules
import re
from collections import defaultdict
from datetime import date
import os

def feedback(dir,conf_dic):
    counts=0
    counts_pair=set()
    n_tfs=0
    tfs_w_eff=0
    feedback_dic={}
    n_tfs = len(os.listdir(os.path.join(dir+"/GUs_secondary_reactions")))
    tfs = [TF for TF in os.listdir(os.path.join(dir+"/GUs_secondary_reactions"))]
    tfs_w_eff = len(set(tfs) & set(conf_dic.keys()))
    tfs_w_eff_set = set(tfs) & set(conf_dic.keys())
    # tfs_wno_eff_set = set(conf_dic.keys()).difference(tfs_w_eff_set)
    tfs_w_feedback = set()
    # print("lens:", tfs_w_eff)
    counts_pair=set()
    for TF in os.listdir(os.path.join(dir+"/GUs_secondary_reactions")):
        feedback_dic[TF] = 0
        if not TF in conf_dic.keys():
            continue
        with open(os.path.join(dir+"/GUs_secondary_reactions/"+TF+"/reactants_products.txt")) as rp:
            for line in rp:
                obj = line.strip().split("\t")[1]
                obj = re.sub('_Ext','',obj)
                if obj in conf_dic[TF]:
                    feedback_dic[TF] += 1
                    counts += 1
                    tfs_w_feedback.add(TF)
                    print(TF,obj,"match")
                    break        
            for line in rp:
                obj = line.strip().split("\t")[1]
                if obj in conf_dic[TF]:
                    counts_pair.add(obj)
    tfs_no_feedback=tfs_w_eff_set.difference(tfs_w_feedback)
    counts=len(tfs_w_feedback)
    return [n_tfs,counts,len(counts_pair),tfs_w_eff,tfs_w_eff_set,tfs_w_feedback,tfs_no_feedback]

# conformationf = "/home/gfemer/GUs/gu_library/TFs_metabolites.txt"
conformationf = "/home/gfemer/GUs/gu_library/TF_metabolite_wgenes.u.txt"
gu_compound_dir = "/home/gfemer/GUs/dictionary/gu_compound_dictionary.KEGG_edited_v2.txt"
gradient_dir = "/home/gfemer/GUs/GUs_geneInteractions_14-05-21_GRADIENT/"

Regulon_dir = "/home/gfemer/GUs/GUs_Regulon10.7_40821"

compound_dic = defaultdict(list)
with open(gu_compound_dir) as comp:
    for line in comp:
        if(re.match("^#", line)):
            continue
        line = line.strip("\n").strip("\r").split("\t")
        compid = line[0]
        compid2 = re.sub("\|","",compid)
        syns = line[2].strip("//").split("//")
        if line[1] not in syns:
            syns.append(line[1])
        syns2 = []
        for syn in syns:
            syn2 = re.sub("\|","",syn)
            if syn2 != syn:
                syns2.append(syn2)
        compound_dic[compid] = compound_dic[compid] + syns + [compid2] + syns2

conformation_dic = defaultdict(list)
n_pairs = 0
with open(conformationf) as conf:
    for line in conf:
        if(re.match("^#", line)):
            continue
        n_pairs += 1 # Each line represents a tf-metabolite pair
        line = line.strip("\n").strip("\r").split("\t")
        # print(line)
        tf = line[0]
        effector = line[1]
        # effector = line[3]
        # effectorKEGG = line[4]
        # conformation_dic[tf] = conformation_dic[tf] + sum([[key] + compound_dic[key] for key, value in compound_dic.items() if any(x in value for x in [effector,effectorKEGG])],[])
        conformation_dic[tf] = conformation_dic[tf] + sum([[key] + compound_dic[key] for key, value in compound_dic.items() if effector in value+[key]],[])
        if(len(sum([[key] + compound_dic[key] for key, value in compound_dic.items() if effector in value+[key]],[]))==0):
            conformation_dic[tf].append(effector)
            

total_tfs_w_eff = len(conformation_dic.keys())
conformation_dic = dict(conformation_dic)
compound_dic = dict(compound_dic)

# table = open(os.path.join(gradient_dir+"GI-140521_Gradient_feedback.txt"),'w')
# table.write("# From search_feedback.py\n")
# table.write("# Date: "+str(date.today())+"\n")
# table.write("# Network gradient directory: "+str(gradient_dir)+"\n")
# table.write("# Network ID\tMin Fe\tMax Fe\tmin -log10(P)\tType of Interactions\tFeedback Counts\tTotal Feedback\tnum TFs with eff\tFeedback from TFs with effectors\tCount of pairs\tProportion of pairs (TF-met) with feedback\n")

#print(len(conformation_dic.keys()))
n_tfs,counts,counts_pair,tfs_w_eff,tfs_w_eff_set,tfs_w_feedback,tfs_no_feedback = feedback(Regulon_dir,conformation_dic)
print("Number of TFs with feedback:",counts)
print("TFs with no feedback:", ", ".join(tfs_no_feedback))

for tf in tfs_no_feedback:
    print(tf,"\t",conformation_dic[tf])

# print(len(conformation_dic.keys()))
# table.write(str(Regulon_dir)+"\t-\t-\t-\t-\t"+str(counts)+"\t"+str(counts/n_tfs)+"\t"+str(tfs_w_eff)+"\t"+str(counts/tfs_w_eff)+"\t"+str(counts_pair)+"\t"+str(counts_pair/n_pairs)+"\n")

# for GU_dir in os.listdir(gradient_dir):
#     m = re.match(r'GI-140521_fe(\d+-\d+)_?p(-?\d\.?\d?)_(\w+)', GU_dir)
#     if not m:
#         continue
#     fe_min=m.group(1).split("-")[0]
#     fe_max=m.group(1).split("-")[1]
#     p=m.group(2)
#     interactions=m.group(3)
#     #print(len(conformation_dic.keys()))
#     n_tfs,counts,counts_pair,tfs_w_eff = feedback(os.path.join(gradient_dir+GU_dir),conformation_dic)
#     #print(tfs_w_eff)
#     #print(str(GU_dir)+"\t"+str(fe_min)+"\t"+str(fe_max)+"\t"+str(p)+"\t"+str(interactions)+"\t"+str(counts)+"\t"+str(counts/n_tfs)+"\t"+str(tfs_w_eff)+"\t"+str(counts/tfs_w_eff)+"\t"+str(counts_pair)+"\t"+str(counts_pair/n_pairs)+"\n")
#     table.write(str(GU_dir)+"\t"+str(fe_min)+"\t"+str(fe_max)+"\t"+str(p)+"\t"+str(interactions)+"\t"+str(counts)+"\t"+str(counts/n_tfs)+"\t"+str(tfs_w_eff)+"\t"+str(0 if tfs_w_eff == 0 else counts/tfs_w_eff)+"\t"+str(counts_pair)+"\t"+str(counts_pair/n_pairs)+"\n")

# table.close()
# print("Finished.")
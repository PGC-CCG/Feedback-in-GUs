#!/usr/bin/env python3

import os
import re
import pandas as pd
from itertools import product
from datetime import date

networkf = "/home/gfemer/GUs/Networks/geneInteractions_14-05-21.txt"
gradient_dir = "/home/gfemer/GUs/Networks/geneInteractions_14-05-21_GRADIENT/"
basename = "GI-140521"
os.mkdir(gradient_dir)

network = pd.read_table(networkf,header=0)
network['tf'] = network['tf'].str.capitalize()

# Filters 
min_list = list(range(0,12,2))
max_list = [6000,2000,500,200,100,50,30]
p_val = [-1,1,1.3,1.5,2,2.5,3]
combinations = list(product(min_list,max_list,p_val))

# known known_weak || known_strong
# known_weak known_weak == 1
# known_strong known_strong == 1
# unknown !(known_weak || known_strong)
# all sin filtro

for pars in combinations:
    # min and max values in max_fe
    min_fe = pars[0]
    max_fe = pars[1]
    p = pars[2]
    # All
    f1 = network.query('max_fe >= @min_fe & @max_fe <= @max_fe & pval >= @p')
    filename = str(gradient_dir)+str(basename)+"_fe"+str(min_fe)+"-"+str(max_fe)+"_p"+str(p)+"_all.txt"
    with open(filename, 'a') as f:
        f.write("# From: create_network_gradient.py on {}\n# Input Network: {}\n# Filters: Min_Fe = {}, Max_Fe = {}, -log(P) = {}\n# TF   Gene\n".format(str(date.today()),networkf,min_fe,max_fe,p))
        f1.to_csv(f,sep="\t",columns=['tf','gene'],header = False, index = False)
    
    # Known weak
    f1 = network.query('max_fe >= @min_fe & max_fe <= @max_fe & pval >= @p & known_weak == 1')
    filename = str(gradient_dir)+str(basename)+"_fe"+str(min_fe)+"-"+str(max_fe)+"_p"+str(p)+"_KnownWeak.txt"
    with open(filename, 'a') as f:
        f.write("# From: create_network_gradient.py on {}\n# Input Network: {}\n# Filters: Min_Fe = {}, Max_Fe = {}, -log(P) = {}, Interactions = Known Weak\n# TF   Gene\n".format(str(date.today()),networkf,min_fe,max_fe,p))
        f1.to_csv(f,sep="\t",columns=['tf','gene'],header = False, index = False)
    
    # Known Strong
    f1 = network.query('max_fe >= @min_fe & max_fe <= @max_fe & pval >= @p & known_strong == 1')
    filename = str(gradient_dir)+str(basename)+"_fe"+str(min_fe)+"-"+str(max_fe)+"_p"+str(p)+"_KnownStrong.txt"
    with open(filename, 'a') as f:
        f.write("# From: create_network_gradient.py on {}\n# Input Network: {}\n# Filters: Min_Fe = {}, Max_Fe = {}, -log(P) = {}, Interactions = Known Strong\n# TF   Gene\n".format(str(date.today()),networkf,min_fe,max_fe,p))
        f1.to_csv(f,sep="\t",columns=['tf','gene'],header = False, index = False)
    
    # Unknown sites
    f1 = network.query('max_fe >= @min_fe & max_fe <= @max_fe & pval >= @p & (not (known_weak | known_strong))')
    filename = str(gradient_dir)+str(basename)+"_fe"+str(min_fe)+"-"+str(max_fe)+"_p"+str(p)+"_UnknownInteractions.txt"
    with open(filename, 'a') as f:
        f.write("# From: create_network_gradient.py on {}\n# Input Network: {}\n# Filters: Min_Fe = {}, Max_Fe = {}, -log(P) = {}, Interactions = Unknown\n# TF   Gene\n".format(str(date.today()),networkf,min_fe,max_fe,p))
        f1.to_csv(f,sep="\t",columns=['tf','gene'],header = False, index = False)
    
    # Only Known sites
    f1 = network.query('max_fe >= @min_fe & max_fe <= @max_fe & pval >= @p & (known_weak | known_strong)')
    filename = str(gradient_dir)+str(basename)+"_fe"+str(min_fe)+"-"+str(max_fe)+"_p"+str(p)+"_KnownInteractions.txt"
    with open(filename, 'a') as f:
        f.write("# From: create_network_gradient.py on {}\n# Input Network: {}\n# Filters: Min_Fe = {}, Max_Fe = {}, -log(P) = {}, Interactions = Known \n# TF   Gene\n".format(str(date.today()),networkf,min_fe,max_fe,p))
        f1.to_csv(f,sep="\t",columns=['tf','gene'],header = False, index = False)



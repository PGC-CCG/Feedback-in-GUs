#!/bin/bash

###################################################
## Set operations between 2 networks             ##
## Author: Georgette Femerling                   ##
## Version: v1                                   ##
## Date: 1-03-2020                               ##
## Description: Script that takes the output     ##
## files from the remodel_reactions.py script 	 ##
## from two networks and compares the two sets 	 ##
## of reactions                                  ##
## Input: [1] Modelled reactions dir1            ##
##        [2] Modelled reactions dir2            ##
##        [3] comparison file name               ##
##        [4] output directory                   ##
## Output: One comparison table                  ##
## 3 files per GU: UniqueA, UniqueB and compare  ##
###################################################

# run as: bash ./perl_scripts/set_operations.sh Regulon10.6/datasets/Modelled_GUs/ PhysicalNt_groups2/datasets/Modelled_GUs/ comp_regphy_groups.txt Comparison

# Before running
# Changes the names of files to first Uppercase agaR -> AgaR
# for i in *; do mv "$i" "${i^}"; done

# Take parameters.
source_data_A=$1 #dir of modelled reactions of first group
source_data_B=$2 #dir of modelled reactions of second group
out_data=$3 #output file name
out_dir=$4 #output directory

#Create output directory
mkdir $out_dir

# Identify shared files. - Obtiene los archivos con nombres iguales
intersection_source_data=$(comm -12 <(ls $source_data_A | sort) <(ls $source_data_B | sort))

# Print header.
A=$(echo "$source_data_A"|cut -f1 -d "/")
B=$(echo "$source_data_B"|cut -f1 -d "/")
echo -e "# Source A: $A\n# Source B: $B\n# Date: $(date)" > ${out_dir}/$out_data
echo -e "tf\tunique_A\tunique_B\tintersection_AB\ttotal" >> ${out_dir}/$out_data

# Iterate shared files.
for item in $intersection_source_data; do

	# Get name of item.
	name=$(echo "$item" | cut -f1 -d "_")
		
	# Sort data.
	group_A=$(sort $source_data_A$item | grep -v "#")
	group_B=$(sort $source_data_B$item | grep -v "#")

	#set operations
	unique_A=$(comm -23 <(echo "$group_A") <(echo "$group_B") | wc -l)
	echo -e "Unicas A \n$(comm -23 <(echo "$group_A") <(echo "$group_B"))" > ${out_dir}/${name}.A.txt
	unique_B=$(comm -13 <(echo "$group_A") <(echo "$group_B") | wc -l)
	echo -e "Unicas B \n$(comm -13 <(echo "$group_A") <(echo "$group_B"))" > ${out_dir}/${name}.B.txt
	intersection_AB=$(comm -12 <(echo "$group_A") <(echo "$group_B") | wc -l)
	total=$(echo "$unique_A + $unique_B + $intersection_AB" | bc)
	echo -e "$(comm --output-delimiter=/ <(echo "$group_A") <(echo "$group_B"))" |  \
	awk '{if ($0 ~ /^\/[^\/]/) {print "PhysicalNt\t"substr($0,2);} else if ($0 ~ /^\/\/[^\/]/) {print "shared\t"substr($0,3);} else {print "Regulon\t"$0;}}' \
	> ${out_dir}/${name}.compare.txt	
		
	# Print results.
	echo -e "$name\t$unique_A\t$unique_B\t$intersection_AB\t$total" >>  ${out_dir}/$out_data

done

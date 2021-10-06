#!/bin/bash

################################################################################
######## Main pipeline to generate GUs from group of genes          ############
######## code used to generate GUs: gu_assembly_groups_enz_fromfile.pl  ########
######## Author: Daniela Ledezma-Tejeida                                ########
######## Adapted by: Georgette Femerling                                ########
######## Date: 03-18-2020                                               ########
######## Arguments: [1] regulated genes sep by //                       ########
######## 03-27-2020 Added step 4. modelling template                    ########
################################################################################
#Run pathway-tools on the background as ~/pathway-tools/pathway-tools -api -lisp

Network=$1
network_dir=$2

#GU_assembly on TFs
mkdir ./${network_dir}/GUs_raw_files
echo "GU assembly"
perl ./bin/gu_assembly_groups_enz_fromfile.pl  $Network ./${network_dir}/GUs_raw_files/ ./${network_dir}/gu_library/

#Compare super-reactions to choose appropriate allowed steps 
mkdir ./${network_dir}/tmp1 ./${network_dir}/tmp2 ./${network_dir}/results
echo "Comparing super-reactions"
perl ./bin/secondary_reaction_comparison.pl ./${network_dir}/GUs_raw_files/ ./${network_dir}/tmp1/ ./${network_dir}/tmp2/ ./${network_dir}/tmp_file.txt ./${network_dir}/results/SRs_comparison.txt
rm -r ./${network_dir}/tmp1 ./${network_dir}/tmp2 ./${network_dir}/tmp_file.txt

#Add super_reactions (9 allowed steps)
mkdir ./${network_dir}/GUs_secondary_reactions
echo "Adding super reactions"
perl ./bin/add_secondary_reactions_v4.pl ./${network_dir}/GUs_raw_files/ ./${network_dir}/GUs_secondary_reactions/ ./${network_dir}/results/counts_secondary_reactions.txt 9

#Modeling GUs to new format
echo "Modeling GUs"
mkdir ./${network_dir}/Modelled_GUs
python ./bin/remodel_reactions.py ./${network_dir}/GUs_secondary_reactions/ ./${network_dir}/Modelled_GUs/

#Enzyme counts
echo "Enzyme counts"
perl ./bin/enzyme_count_v2.pl ./${network_dir}/GUs_secondary_reactions/ ./${network_dir}/results/enzyme_counts.txt

#Catalysis reactions counts (for metadata)
echo "Counting catalysis reactions"
perl ./bin/get_catalysis_reactions.pl ./${network_dir}/GUs_secondary_reactions/ ./${network_dir}/results/catalysis_reactions.txt

#Get all short paths (for components)
echo "Getting short paths"
perl ./bin/get_short_paths_v6.pl ./${network_dir}/GUs_secondary_reactions/ ./${network_dir}/results/all_short_paths.txt 0

#Get components
echo "Getting Components"
perl ./bin/find_components.pl ./${network_dir}/results/all_short_paths.txt ./${network_dir}/results/components.txt

#Get short paths (for connectivity)
echo "Getting Short Paths for connectivity"
perl ./bin/get_short_paths_v6.pl ./${network_dir}/GUs_secondary_reactions/ ./${network_dir}/results/short_paths_conn.txt 1

#Get connectivity 
echo "Getting connectivity"
perl ./bin/get_connectivity_v1.pl ./${network_dir}/results/short_paths_conn.txt ./${network_dir}/GUs_secondary_reactions/ ./${network_dir}/results/connectivity.txt

#Get feedback
echo "Getting feedback"
perl ./bin/find_feedback.pl ./${network_dir}/results/all_short_paths.txt ./${network_dir}/GUs_secondary_reactions/ ./${network_dir}/results/feedback.txt

#Get summary
echo "Writting summary"
perl ./bin/get_summary.pl ./${network_dir}/results/ ./${network_dir}/gu_library/Transcription_Factors.txt ./${network_dir}/results/summary.txt


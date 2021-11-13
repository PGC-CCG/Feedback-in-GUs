#!/bin/bash

################################################################################
######## Main pipeline to generate GUs from group of genes          ############
######## code used to generate GUs: gu_assembly_groups_enz_fromfile.pl  ########
######## Author: Daniela Ledezma-Tejeida                                ########
######## Adapted by: Georgette Femerling                                ########
######## Date: 03-18-2020                                               ########
######## Arguments: [1] regulated genes sep by //                       ########
########            [2] directory to write results                      ########
######## 03-27-2020 Added step 4. modelling template                    ########
################################################################################

Netk=$1
output_dir=${2%/} #Remove the last slash if there was

echo -e "............... Loading Pathway tools ..................\n"
screen -d -m -S Ptools bash -c 'pathway-tools -api -lisp'
sleep 5


mkdir $output_dir
mkdir ${output_dir}/tmp

# Create dataset
echo -e ".......... Creating Dataset ................. \n" 
name=$(basename $Netk)
perl bin/create_datasets.pl $Netk $name ${output_dir}/tmp/
Network=${output_dir}/tmp/${name}"_general_regulons.txt"

#GU_assembly on TFs
mkdir ./${output_dir}/GUs_raw_files
echo -e "................ GU assembly ........................... \n"
perl ./bin/gu_assembly_groups_enz_fromfile.pl $Network ./${output_dir}/GUs_raw_files/ ./gu_library/

#Compare super-reactions to choose appropriate allowed steps 
mkdir ./${output_dir}/tmp1 ./${output_dir}/tmp2 ./${output_dir}/results
echo -e "................ Comparing super-reactions ....................... \n"
perl ./bin/secondary_reaction_comparison.pl ./${output_dir}/GUs_raw_files/ ./${output_dir}/tmp1/ ./${output_dir}/tmp2/ ./${output_dir}/tmp_file.txt ./${output_dir}/results/SRs_comparison.txt
rm -r ./${output_dir}/tmp1 ./${output_dir}/tmp2 ./${output_dir}/tmp_file.txt

#Add super_reactions (9 allowed steps)
mkdir ./${output_dir}/GUs_secondary_reactions
echo -e ".....................Adding super reactions ...................... \n"
perl ./bin/add_super_reactions_v5.pl ./${output_dir}/GUs_raw_files/ ./gu_library/ ./${output_dir}/GUs_secondary_reactions/ ./${output_dir}/results/counts_secondary_reactions.txt 9

#Modeling GUs to linear format
echo -e "..................... Modeling GUs ............................. \n"
mkdir ./${output_dir}/Modelled_GUs
python ./bin/remodel_reactions.py ./${output_dir}/GUs_secondary_reactions/ ./${output_dir}/Modelled_GUs/

#Getting connectivity
echo -e "............................ Calculating Connectivity ................... \n"
perl ./bin/connectivity_from_gufiles_permissive_enx_v1.pl ./${output_dir}/GUs_secondary_reactions/ ./${output_dir}/results/connectivity.txt

#Getting Pathways and DeltaG from reactions
echo -e "....................... Getting pathways and DeltaG in reactions ..................... \n"
perl ./bin/get_pathways_deltaG_from_rxs.pl ./${output_dir}/GUs_secondary_reactions/ ./gu_library/ ./${output_dir}/results/

#Add objects ID
echo -e ".......................Adding IDs to objects in GU........................................\n"
mkdir ./${output_dir}/GUs_secondary_reactions_IDs
perl ./bin/add_objects_ID_v1.pl ./${output_dir}/GUs_secondary_reactions ./${output_dir}/GUs_secondary_reactions_IDs

echo -e "........................ GU Assembly done ................................\n"
rm -r ${output_dir}/tmp/
screen -X -S Ptools quit
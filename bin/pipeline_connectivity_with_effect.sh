####
#### Get connectivity considering regulatory effect on enzymes
####


## Get short paths considering effect. ##
# Raw GUs are necessary to avoid secondary reactions to be counted as individual components. It does not happen in common short paths because they are (by definition) connected to a reaction, but here it might not be the case and it alters connectivity calculation.
# Short paths for connectivity must have a limit of 1 (to only include reactions actually connected). Components need a limit of 0 to consider components of only 1 reaction.
# Secondary Reactions are useful for connectivity to link individual reactions.
echo -e "Getting short paths"
perl ../perl_scripts/get_short_paths_v5.pl ./datasets/GUs_secondary_reactions/ ./datasets/connectivity_with_effect/short_paths_activation.txt 1 +
perl ../perl_scripts/get_short_paths_v5.pl ./datasets/GUs_raw_files/ ./datasets/connectivity_with_effect/paths_activation_components.txt 0 +
perl ../perl_scripts/get_short_paths_v5.pl ./datasets/GUs_secondary_reactions/ ./datasets/connectivity_with_effect/short_paths_inhibition.txt 1 -
perl ../perl_scripts/get_short_paths_v5.pl ./datasets/GUs_raw_files/ ./datasets/connectivity_with_effect/paths_inhibition_components.txt 0 -
 
# Enzyme count
echo -e "Counting Enzymes"
perl ../perl_scripts/enzyme_count_effect.pl ../datasets/GUs_secondary_reactions/ ./datasets/connectivity_with_effect/enzyme_counts.txt

# Get connected enzymes
echo -e "Getting connected enzymes"
perl ../perl_scripts/get_connectivity_weffect.pl ./datasets/connectivity_with_effect/short_paths_activation.txt ./datasets/GUs_secondary_reactions ./datasets/connectivity_with_effect/enzyme_counts.txt ./datasets/connectivity_with_effect/connectivity_activation.txt +
perl ../get_connectivity_weffect.pl ./datasets/connectivity_with_effect/short_paths_inhibition.txt ./datasets/GUs_secondary_reactions/ ./datasets/connectivity_with_effect/enzyme_counts.txt ./datasets/connectivity_with_effect/connectivity_inhibition.txt -

#Get components
echo -e "Getting components"
perl ../perl_scripts/find_components.pl ./datasets/connectivity_with_effect/paths_activation_components.txt ./datasets/connectivity_with_effect/activation_components.txt
perl ../perl_scripts/find_components.pl ./datasets/connectivity_with_effect/paths_inhibition_components.txt ./datasets/connectivity_with_effect/inhibition_components.txt 

#Connectivity calculation
echo -e "Connectivity"
perl ../perl_scripts/connectivity_calculator_weffect.pl ./datasets/connectivity_with_effect/ ./datasets/connectivity_with_effect/connectivity_summary.txt

#Correr pathway-tools
# El ReRun de Dani fue con 189 TFs
#~/pathway-tools/pathway-tools -api -lisp
#system(pathway-tools -api -lisp)

#GU_assembly on TFs
mkdir ./datasets/GUs_raw_files
echo "GU assembly"
perl ../perl_scripts/gu_assembly_v5_GF.pl d ./datasets/gu_library/Transcription_Factors.txt ./datasets/GUs_raw_files/ ./datasets/gu_library/

#Compare super-reactions to choose appropriate allowed steps 
mkdir ./datasets/tmp1 ./datasets/tmp2 ./datasets/results
echo "Comparing super-reactions"
perl ../perl_scripts/secondary_reaction_comparison.pl ./datasets/GUs_raw_files/ ./datasets/tmp1/ ./datasets/tmp2/ ./datasets/tmp_file.txt ./datasets/results/SRs_comparison.txt
rm -r ./datasets/tmp1 ./datasets/tmp2 ./datasets/tmp_file.txt

#Add super_reactions (9 allowed steps)
mkdir ./datasets/GUs_secondary_reactions
echo "Adding super reactions"
perl ../perl_scripts/add_secondary_reactions_v4.pl ./datasets/GUs_raw_files/ ./datasets/GUs_secondary_reactions/ ./datasets/results/counts_secondary_reactions.txt 9

#Enzyme counts
echo "Enzyme counts"
perl ../perl_scripts/enzyme_count_v2.pl ./datasets/GUs_secondary_reactions/ ./datasets/results/enzyme_counts.txt

#Catalysis reactions counts (for metadata)
echo "Counting catalysis reactions"
perl ../perl_scripts/get_catalysis_reactions.pl ./datasets/GUs_secondary_reactions/ ./datasets/results/catalysis_reactions.txt

#Get all short paths (for components)
echo "Getting short paths"
perl ../perl_scripts/get_short_paths_v6.pl ./datasets/GUs_secondary_reactions/ ./datasets/results/all_short_paths.txt 0

#Get components
echo "Getting Components"
perl ../perl_scripts/find_components.pl ./datasets/results/all_short_paths.txt ./datasets/results/components.txt

#Get short paths (for connectivity)
echo "Getting Short Paths for connectivity"
perl ../perl_scripts/get_short_paths_v6.pl ./datasets/GUs_secondary_reactions/ ./datasets/results/short_paths_conn.txt 1

#Get connectivity 
echo "Getting connectivity"
perl ../perl_scripts/get_connectivity_v1.pl ./datasets/results/short_paths_conn.txt ./datasets/GUs_secondary_reactions/ ./datasets/results/connectivity.txt

#Get feedback
echo "Getting feedback"
perl ../perl_scripts/find_feedback.pl ./datasets/results/all_short_paths.txt ./datasets/GUs_secondary_reactions/ ./datasets/results/feedback.txt

#Get summary
echo "Writting summary"
perl ../perl_scripts/get_summary.pl ./datasets/results/ ./datasets/gu_library/Transcription_Factors.txt ./datasets/results/summary.txt


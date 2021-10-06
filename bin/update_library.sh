mkdir gu_library/
cd gu_library/

echo -e "............... Loading Pathway tools ..................\n"
screen -d -m -S Ptools bash -c 'pathway-tools -api -lisp'
sleep 5

echo -e "................. Updating gene links ..............\n"
wget https://ecocyc.org/gene-links.dat
mv gene-links.dat gene_links.txt

echo -e "........... Updating Compound dictionary ............\n"
perl ../bin/create_compound_dictionary.pl ./

echo -e "............ Updating gene-reactions flat file ............\n"
perl ../bin/gene_rxns_flatfiles.pl gene_links.txt gu_compound_dictionary.txt rxn_rpd_links.txt

echo -e "............. Updating pathway-reaction links file ...........\n"
perl ../bin/create_pwy_rxn_links.pl ./

echo -e "............. Downloading TF Conformation from Regulon DB ........... \n"

cd ..
echo -e "\n\n............ GU ibrary files Updated .................. \n\n"
echo -e "........ Done ........ \n"
screen -X -S Ptools quit
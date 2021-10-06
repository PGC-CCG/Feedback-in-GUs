#################################################################################
###### 06-25-2020                                                          ######
###### Re-written by: Georgette Femerling                                  ######
###### Create pwd_rxn_links.txt                                            ######
###### Creates input gu_library file that contains all Ecoli pathways and  ######
###### their reactions                                                     ######
###### Argument usage:                                                     ######
###### [0] - path to gu_library directory                                  ######
###### Output:                                                             ######
###### pwy_rxn_links.txt                                                   ######
###### [1] - pathway ID                                                    ######
###### [2] - reaction ID                                                   ######
###### History:                                                            ######
#################################################################################

#system("pathway-tools -api -lisp");

#Load modules
use perlcyc;
$cyc = perlcyc -> new("ECOLI");

$time=localtime();

#Get input arguments
$dir_out=$ARGV[0];
chomp($dir_out);

#Open outfiles
if($dir_out){
    $pwyrxn=$dir_out . "pwy_rxn_links.txt";
    open(OUT,">$pwyrxn") || die "Cannot open OUT file at $pwyrxn.\n";
}
print OUT"# From $0 on $time\n";

@pathways = $cyc -> all_pathways();
foreach $path (@pathways){
     @rxpath = $cyc -> get_reaction_list($path);
     foreach $rx(@rxpath){
          print OUT"$path\t$rx\n"
     }
}

close(OUT);

print "................Done..................\n";
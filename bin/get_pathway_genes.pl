#################################################################################
###### 11-08-2020                                                          ######
###### Author: Georgette Femerling                                         ######
######                                                                     ######
###### Gets genes of pathways                                              ######
###### Argument usage:                                                     ######
###### [0] - path to gu_library directory                                  ######
###### Output:                                                             ######
######  1 - pathway_genes.txt :                                            ######
######      [1] Pathway                                                    ######
######      [2] genes                                                      ######
#################################################################################

#Load modules
use perlcyc;
$cyc = perlcyc -> new("ECOLI");

#get arguments
$time=localtime();
$dir_out=$ARGV[0];
chomp($dir_out);

$gene_links=$dir_out . "gene_links.txt";

#create output files
$output=$dir_out . "pathway_genes.txt";
$output2=$dir_out . "pathway_names.txt";

#Open outfiles
open(OUT,">$output") || die "Cannot open OUT file at $output.\n";
open(OUT2,">$output2") || die "Cannot open OUT2 file at $output2.\n";

#Print output files headers
print OUT"# From $0 on $time\n# Pathway\tGenes\n";
print OUT2"# From $0 on $time\n# Pathway ID\t Common Name\n";

#Get all ecoli pathways 
@pathways = $cyc -> all_pathways();

my %pwys;

foreach $pwy (@pathways){
    #print OUT"$pwy\t";
    $pw_name = $cyc -> get_slot_value($pwy, "COMMON-NAME");
    @genes = $cyc -> genes_of_pathway($pwy);
    foreach $id(@genes){
        $gname=&IDS($id);
        chomp($gname);
        push (@{$pwys{$pwy}}, $gname);
    }
    print OUT "$pwy\t";
    print OUT join(",", @{$pwys{$pwy}}),"\n";
    print OUT2 "$pwy\t$pw_name\n";
}

close(OUT);
close(OUT2);

####################

sub IDS {
 
  local($id,$gene)=($_[0],"");
  open(IDS, $gene_links) || die "Cannot open IDS at $gene_links.\n";
  while(<IDS>){
    if($_=~/^$id\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t([^\t]+)$/i){
      $gene=$1;
      close(IDS);
      return($gene);
    }
  }
  close(IDS);
  print"WARNING 2! Gene $id is not present in gene_links.\n";
  return($id);
}

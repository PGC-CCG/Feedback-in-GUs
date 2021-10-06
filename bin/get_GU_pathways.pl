#################################################################################
###### 21-06-16                                                            ######
###### Get pathways of GU                                                  ######
###### Finds pathways per gene per GU                                      ######
###### Argument usage:                                                     ######
###### [0] - /datasets/GOs_Pathways_in_GENSORUnits/genes_in_GUs.txt        ######
###### [1] - /datasets/input_files/gene_name_links.txt file                ######
###### [2] - output file                                                   ######
###### Output:                                                             ######
###### [0] > GU name                                                       ######
###### [1] > Total pathways                                                ######
###### [2] > Regulon size                                                  ######
###### [3] > Genes without pathways                                        ######
###### [4] > All uniq pathways (number of times repeated)                  ######
#################################################################################

$time=localtime();
use perlcyc;
$cyc = perlcyc -> new("ECOLI");

#Get args
$reg_file=$ARGV[0];
$gene_links=$ARGV[1];
$outfile=$ARGV[2];
chomp($outfile);

#Print OUT header
open(OUT,">$outfile") || die "Cannot open OUT at $outfile.\n";
print OUT"#From $0 on $time using: #  $reg_file\n# [1] > GU name\n# [2] > Total pathways\n# [3] > Regulon size\n# [4] > Genes without pathway\n# [5] > Individual pathways (number of ocurrences)\nGU name\tTotal pathways\tRegulon size\tGenes without pathway\tIndividual pathways (number of ocurrences)\n";

#Open IN file
open(IN,$reg_file) || die "Cannot open IN file at $reg_file.\n";
while(<IN>){
  if($_=~/^([^\t]+)\t([^\t]+)$/){
    $TF=$1;
    $genes=$2;
    @genes=split(/\/\//,$genes);
    undef %pwys;
    @pwys=();
    $no_pwy=0;
    foreach $g (@genes){
      $gid=&IDS($g);
      @pwys=$cyc->pathways_of_gene($gid);
      foreach $pw (@pwys){
	$pwys{$pw}++;
      }
      if(!(@pwys)){
	$no_pwy++;
      }
    }
    @kys=keys %pwys;
    $sizek=@kys;
    $sizeg=@genes;
    print OUT"$TF\t$sizek\t$sizeg\t$no_pwy\t";
    foreach $p (keys %pwys){
      print OUT"$p($pwys{$p})//";
    }
    print OUT"\n";
  }
}
close(IN);
close(OUT);



###########################
### IDS Function
###
### Turn gene names into Ecocyc IDS
### -Input arguments: gene name
### -Returned values: Ecocyc ID
###########################



sub IDS {
 
  local($name,$id)=($_[0],"");
  open(IDS, $gene_links) || die "Cannot open IDS at $gene_links.\n";
  while(<IDS>){
    if($_=~/^([^\t]+)\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t$name$/i){
      $id=$1;
      close(IDS);
      return($id);
    }
  }
  close(IDS);
  print"WARNING 2! Gene $name is not present in gene_links.\n";
  return($name);
}








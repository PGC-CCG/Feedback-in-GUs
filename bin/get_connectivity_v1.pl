#################################################################################
###### 26-11-15                                                            ######
###### Calculate connectivity                                              ######
###### Calculates connectivity from get_short_paths_vX.pl                  ######
###### Arguments:                                                          ######
###### [0]> get_short_paths_vX.pl output                                   ######
###### [1]> GU_assembly output dir                                         ######
###### [2]> Output file                                                    ######
###### Output:                                                             ######
###### [1]> GU name                                                        ######
###### [2]> Total enzymes with connectivity                                ######
###### [3]> Enzymes with connectivity                                      ######
###### History:                                                            ######
###### 3/04/16 > V1 > Complexes are quantified as a single enzyme          ######
#################################################################################

$time=localtime();

#Get args
$infile=$ARGV[0];
$gu_dir=$ARGV[1];
$outfile=$ARGV[2];
chomp($outfile);

undef %connectivity;

#Open files
open(IN,$infile) || die "Cannot open IN file at $infile.\n";
open(OUT,">$outfile") || die "Cannot open OUT file at $outfile\n";
print OUT"# From get_connectivity.pl on $time\n# Input = $infile / $gu_dir \n";

while(<IN>){
  if($_=~/^([^\t]+)\t([^\t]+)\t/){   #read short_paths
    $TF=$1;                          
    $path=$2;
    
    #Open GU_assembly files
    $mods=$gu_dir . $TF . "/modification.txt";

    @rxns=split(/\/\//,$path);           #split TF reactions

    foreach $rx (@rxns){

      #Find modifiers
      open(MOD,$mods) || die "Cannot open MOD at $mods.\n";    #look for catalysts
      while(<MOD>){
	if($_=~/^$rx\tCATALYSIS\t([^\t]+)$/i){
	  $mod=$1;
	  chomp($mod);
	  #Monomers
	  #Evaluate if modifier is already in vector
	  $rp=&repeats($mod,\@{$connectivity{$TF}});
	  if(!($rp)){
	    push(@{$connectivity{$TF}},$mod);
	  }
	}
      }
      close(MOD);

    }
    ##Print TFs with no connectivity for easier data handling of all GUs
  }elsif($_=~/^([^\t]+)\t\t\t/){
    print OUT"$1\t0\t\n";
  }
}
close(IN);


#Print total modifiers.

foreach $ft (keys %connectivity){
  $size=(@{$connectivity{$ft}});
  print OUT"$ft\t$size\t";
  foreach $m (@{$connectivity{$ft}}){
    print OUT"$m,";
  }
  print OUT"\n";
}




###########################
### repeats Function
###
### Finds if a value is already present in a vector
### -Input arguments:
###  [0] - value
###  [1] - vector
### -Returned values: 0 if absent; 1 if present
###########################


sub repeats{

  local($val,@vector,$elem)=($_[0],@{$_[1]},"");
  
  foreach $elem (@vector){
    if($elem eq $val){
      return(1);
    }
  }

  return(0);
}






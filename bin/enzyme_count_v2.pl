#################################################################################
###### 04-09-14                                                            ######
###### Count enzymes per GU                                                ######
###### Counts number of proteins that have a CATALYSIS reaction associated ######
###### in modifications.txt                                                ######
###### Arguments:                                                          ######
###### [0]> GU_assembly output folder path                                 ######
###### [1]> Output file                                                    ######
###### Output:                                                             ######
###### [1]> GU name                                                        ######
###### [2]> Number of enzymes                                              ######
###### [3]> Protein names separated by commas                              ######
###### History:                                                            ######
###### v1 > 17/11/15 > Code rewritten                                      ######
###### v2 > 03/04/16 > Complexes are quantified as single enzymes          ######
#################################################################################

$time=localtime();

## Get arguments
$GU_path=$ARGV[0];
$outfile=$ARGV[1];
chomp($outfile);

##Open outfile
open(OUT,">$outfile") || die "Cannot open OUT at $outfile.\n";
print OUT"# From enzyme_count_v1.pl on $time\n# Input = $GU_path\n";

## Get GU dirs
opendir(DIR,$GU_path) || die "Cannot open $GU_path.\n";
@GU_files=readdir(DIR);
foreach $file (@GU_files){
  if(!($file=~/^\./)){

    #Get modifications file
    $modifications=$GU_path . $file . "/modification.txt";

    #Initialize variables
    @enzymes=();

    #Find CATALYSIS modifications
    open(MOD,$modifications) || die "Cannot open MOD at $modifications.\n"; 
    while(<MOD>){
      if($_=~/^re\d+\tCATALYSIS\t([^\t]+)$/i){
	$prot=$1;
	chomp($prot);
	$flag=0;
	foreach $rp (@enzymes){
	  if($prot eq $rp){
	    $flag=1;
	    last;
	  }
	}
	if(!($flag)){
	  push(@enzymes,$prot);
	}	  
      }
    }

    $size=@enzymes;
    print OUT"$file\t$size\t";
    foreach $ez (@enzymes){
      print OUT"$ez,";
    }
    print OUT"\n";
  }
}
closedir(DIR);

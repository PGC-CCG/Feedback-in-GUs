#################################################################################
###### 06-21-2020                                                          ######
###### Author: Georgette Femerling                                         ######
######                                                                     ######
###### Extracts the DeltaG as well as the pathways of a reaction in a GU   ######
###### Per GU extracts the pathways present and calculates the total number######
###### Extracts pathways per GU considering only reactions with a deltaG   ######
######  less than 0                                                        ######
###### Argument usage:                                                     ######
###### [0] - path to GUs_secondary_reactions directory                     ######
###### [1] - path to gu_library directory                                  ######
###### [2] - path to output directory (results dir)                        ######
###### Output:                                                             ######
######  1 - rx_pathways_deltaG.txt :                                       ######
######      [1] GU                                                         ######
######      [2] GU rx number                                               ######
######      [3] rx ID                                                      ######
######      [4] pathways                                                   ######
######      [5] number of pathways                                         ######
######      [6] DeltaG                                                     ######
######      [7] direction                                                  ######
######  2 - gu_pathways.txt                                                ######
######      [1] GU                                                         ######
######      [2] Ecocyc pathways in GU                                      ######
######      [3] Number of pathways                                         ######
######  3 - gu_pathways_deltaG0.txt                                        ######
######      [1] GU                                                         ######
######      [2] Ecocyc pathways in GU, in rxs with deltaG < 0              ######
######      [3] Number of pathways                                         ######
###### History:                                                            ######
#################################################################################

#Load modules
use perlcyc;
$cyc = perlcyc -> new("ECOLI");

#get arguments
$time=localtime();
$in_dir=$ARGV[0];
$gulib=$ARGV[1];
$dir_out=$ARGV[2];
chomp($dir_out);

#Open pwy_rx_links
$pwy_link=$gulib . "pwy_rxn_links.txt";

#Create outputfiles
$output=$dir_out . "rx_pathways_deltaG.txt";
$output2=$dir_out . "gu_pathways.txt";
$output3=$dir_out . "gu_pathways_deltaG0.txt";

#Open outfiles
open(OUT,">$output") || die "Cannot open OUT file at $output.\n";
open(OUT2,">$output2") || die "Cannot open OUT2 file at $output2.\n";
open(OUT3,">$output3") || die "Cannot open OUT3 file at $output3.\n";


#Print output files headers
print OUT"# From $0 on $time\n# Indir = $in_dir\nGU\treaction\treaction ID\tpathways\tnumber of pathways\tdeltaG\tdirection\n";
print OUT2"# From $0 on $time\n# Indir = $in_dir\nGU\tEcocyc pathways in GU\tnumber of pathways\n";
print OUT3"# Ecocyc pathways in reactions with negative DeltaG in GUs\n# From $0 on $time\n# Indir = $in_dir\nGU\tEcocyc pathways in GU\tnumber of pathways\n";

#Open IN dir
opendir(DIR,$in_dir) || die "Cannot open DIR at $in_dir.\n";
@GU_files=readdir(DIR);

foreach $file (sort @GU_files){
  if(!($file=~/^\./)){

    print"\n**$file\n";
    $TF=$file;
    $c_rxn=0;
    my %all_pathways;
    my %deltaG0_pathways;

    #Get file names from GU folder
    $o_reactions=$in_dir . $file . "/reactions.txt";

    #Save all reactions
    open(RXN,$o_reactions) || die "Cannot open RXN file at $reactions.\n";
    while(<RXN>){
      @pathways=();
      #Get reaction
      if($_=~/^re(\d+)\t/){
	      $c_rxn=$1;
      }
      if($_=~/^re\d+\t[^\t]+\t[^\t]+\t([^\t]+)$/){   #get ECOCYC reaction ID
        $rx=$1;
        chomp($rx);

        #Get DeltaG from rx
        $DeltaG = $cyc -> get_slot_value($rx, "GIBBS-0");

        #Get rx pathways
        open(LNK,$pwy_link) || die "Cannot open LNK at $pwy_link.\n";
          while(<LNK>){
               if($_=~/^([^\t]+)\t$rx$/){   
	              push(@pathways,$1);      #get pathways of reaction
                foreach $pt(@pathways){
                    if (!exists($all_pathways{$pt})){
                      $all_pathways{$pt} = 1;
                    }
                    if (($DeltaG < 0) && (!exists($deltaG0_pathways{$pt}))){
                      $deltaG0_pathways{$pt} = 1;
                    }
               }
	          }
          }
        close(LNK);
    
        
        #Print first output
        if(@pathways > 0){
          print OUT"$TF\tre$c_rxn\t$rx\t";
          print OUT join(",",@pathways),"\t";
          print OUT scalar(@pathways),"\t";
          print OUT"$DeltaG\t";

          #Get Rx direction
          @f_phyr = $cyc -> get_slot_values($rx,"REACTION-DIRECTION");
          print OUT"$f_phyr[0]\n";
        }
      }
      
    }
    close(RXN);

    #print second output
    if (keys %all_pathways > 0){
      print OUT2 "$TF\t";
      print OUT2 join(",", keys %all_pathways),"\t";
      print OUT2 scalar(keys %all_pathways),"\n";
    }

    #print third output
    if (keys %deltaG0_pathways > 0){
      print OUT3 "$TF\t";
      print OUT3 join(",", keys %deltaG0_pathways),"\t";
      print OUT3 scalar(keys %deltaG0_pathways),"\n";
    }
  } 
}
close(OUT);
close(OUT2);
close(OUT3);

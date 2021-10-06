#################################################################################
###### 07-12-15                                                            ######
###### Secondary-reactions comparison                                      ######
###### Runs add_secondary_reactions_v4.pl allowing 0-25 intermediate metabolites. ###
###### Returns matrix with added super-reactions                           ######
###### Argument usage:                                                     ######
###### [0] - GU_assembly raw data dir                                      ######
###### [1] - Temporary IN dir                                              ######
###### [2] - Temporary OUT dir                                             ######
###### [3] - Temporary OUT file                                            ######
###### [4] - Matrix OUT file                                               ######
###### History:                                                            ######
#################################################################################

$time=localtime();

#Get args
$gu_dir=$ARGV[0];
$temp_in=$ARGV[1];
$temp_out=$ARGV[2];
$temp_file=$ARGV[3];
$outfile=$ARGV[4];
chomp($outfile);

#Print OUT file header
open(OUT,">$outfile") || die "Cannot open OUT file at $outfile.\n";
print OUT"# From super_reaction_comparison.pl on $time\n";

#Open DIR
opendir(DIR,$gu_dir) || die "Cannot open DIR at $gu_dir.\n";
@GU_files=readdir(DIR);
foreach $file (@GU_files){
  if(!($file=~/^\./)){

    print"\n**$file\n";
    print OUT"$file\t";

    $infile=$gu_dir . $file;
    $temp_in_dir=$temp_in . $file;

    #Copy file into IN dir
    system("cp -R $infile $temp_in");
 ###Call function for all metabolite limits
    for($cont=1;$cont<=25;$cont++){
      system("perl ./bin/add_secondary_reactions_v4.pl $temp_in $temp_out $temp_file $cont");
      open(TMP,$temp_file) || die "Cannot open TMP at $temp_file.\n";
      while(<TMP>){
	if((!($_=~/^#/)) && ($_=~/^([^\t]+)\t(\d+)/)){
	  $tf=$1;
	  $srs=$2;
	  chomp($srs);
	  if($tf eq $file){
	    print OUT"$srs\t";
	  }else{
	    print OUT"$tf\t";
	  }
	}
      }
      close(TMP);
    }
    system("rm -rf $temp_in_dir");
    print OUT"\n";
   
  }
}
closedir(DIR);

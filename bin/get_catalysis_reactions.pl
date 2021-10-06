#################################################################################
###### 19-04-16                                                            ######
###### Get Catalysis Reactions                                             ######
###### Creates a file with total catalysis reactions per GU (omitting SRs) ######
###### [0] - TF                                                            ######
###### [1] - Total catalysis reactions                                     ######
###### Argument usage:                                                     ######
###### [0] - GU_assembly directory                                         ######
###### [1] - Output file                                                   ######
#################################################################################

$time=localtime();

#Get args
$in_dir=$ARGV[0];
$outfile=$ARGV[1];
chomp($outfile);

#Open & print header of outfile
open(OUT,">$outfile") || die "Cannot open OUT file at $outfile.\n";
print OUT"# From get_catalysis_reactions.pl on $time\n";

#Open IN dir
opendir(DIR,$in_dir) || die "Cannot open DIR at $in_dir.\n";
@GU_files=readdir(DIR);
foreach $file (@GU_files){
  if(!($file=~/^\./)){

    print"\n**$file\n";
    $TF=$file;

    #Get reactions.txt file path
    $reactions=$in_dir . $TF . "/reactions.txt";

    #Reset counter
    $cont=0;

    #Open reactions files
    open(RXN,$reactions) || die "Cannot open RXN file at $reactions\n";
    while(<RXN>){
      if(($_=~/^re\d+\tSTATE_TRANSITION\t[^\t]+\t/) || ($_=~/^re\d+\tTRANSPORT\t[^\t]+\t/)){
	$cont++;
      }
    }
    close(RXN);

    #Print total reactions
    print OUT"$TF\t$cont\n";
  }
}
closedir(DIR);

#################################################################################
###### 17-11-15                                                            ######
###### Find Feedback                                                       ######
###### Finds feedback in get_short_paths_vX.pl output                      ######
###### [0] - get_short_paths_vX.pl output file                             ######
###### [1] - gu_assembly_vX.pl output dir (to get effectors)               ######
###### [2] - Output file path                                              ######
###### Improvements over get_feedback.pl:                                  ######
###### - Uses get_short_path.pl output instead of re-calculating from      ######
######   gu_assembly.pl                                                    ######
###### - Considers proteins as effectors                                   ######
###### History:                                                            ######
#################################################################################

$time=localtime();

#Get args
$infile=$ARGV[0];
$indir=$ARGV[1];
$outfile=$ARGV[2];
chomp($outfile);

#Open outfile
open(OUT,">$outfile") || die "Cannot open OUT at $outfile.\n";
print OUT"# From find_feedback.pl on $time\n# Input = $infile / $indir \n";

#Open IN file
open(IN,$infile) || die "Cannot open IN file at $infile\n";
while(<IN>){
  if($_=~/^([^\t]+)\t[^\t]+\t([^\t]+)\t/){
    $TF=$1;
    $path=$2;
    
    ###Find effector
    @effectors=();
    #Open complex file
    $complex=$indir . $TF . "/complexes.txt";
    open(CPX,$complex) || die "Cannot open CPX at $complex.\n";
    while(<CPX>){
      if($_=~/^csa\d+\t$TF-[^\t]+\t([^\t]+)/i){
	$eff=$1;
	chomp($eff);
	if($eff ne $TF){
	  push(@effectors,$eff);
	}
      }
    }
    close(CPX);

    #Search all effectors in path
    foreach $ef (@effectors){
      $q_ef=quotemeta($ef);
      if($path=~/$q_ef/i){
	 print OUT"$TF\t$ef\t$path\n";
       }
    }
  }
}
close(IN);

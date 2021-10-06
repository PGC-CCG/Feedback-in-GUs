#################################################################################
###### 19-04-16                                                            ######
###### Obtain Summary                                                      ######
###### Creates a file with columns for data analysis                       ######
###### [0] - GU                                                            ######
###### [1] - Total enzymes                                                 ######
###### [2] - Enzymes with connectivity                                     ######
###### [3] - Components                                                    ######
###### [4] - Connectivity                                                  ######
###### [5] - Total Catalysis reactions                                     ######
###### [6] - Secondary Reactions                                           ######
###### [7] - Feedback                                                      ######
###### Argument usage:                                                     ######
###### [0] - Directory with results files                                  ######
###### [1] - File with TFs to get metadata                                 ######
###### [2] - Output file                                                   ######
#################################################################################

$time=localtime();

#Get args
$in_dir=$ARGV[0];
$infile=$ARGV[1];
$outfile=$ARGV[2];
chomp($outfile);

#Get files in in_dir
$components=$in_dir . "components.txt";
$connectivity=$in_dir . "connectivity.txt";
$feedback=$in_dir . "feedback.txt";
$srs=$in_dir . "counts_secondary_reactions.txt";
$enzyme=$in_dir . "enzyme_counts.txt";
$reactions=$in_dir . "catalysis_reactions.txt";

#Open OUT file and print header
open(OUT,">$outfile") || die "Cannot open OUT at $outfile.\n";
print OUT"#From get_metadata.pl on $time using:\n# $in_dir\n# $infile\nGU\tTotal Enzymes\tEnzymes with connectivity\tComponents\tConnectivity\tTotal Catalysis reactions\tSecondary reactions\tFeedback\n";

##Get TFs
open(IN,$infile) || die "Cannot open IN file at $infile.\n";
while(<IN>){

  if($_=~/^([^\t]+)$/){
    $TF1=$1;
    chomp($TF1);
    $TF=quotemeta($TF1);


    undef %data;

    #Get total enzymes
    $slot=1;
    open(TEN,$enzyme) || die "Cannot open TEN at $enzyme.\n";
    while(<TEN>){
      if($_=~/^$TF\t(\d+)\t/){
	$data{$slot}=$1;
	last;
      }
    }
    close(TEN);

    #Get enzymes with connectivity
    $slot=2;
    open(CONN,$connectivity) || die "Cannot open $connectivity\n";
    while(<CONN>){
      if($_=~/^$TF\t(\d+)\t/){
	$data{$slot}=$1;
	last;
      }
    }
    close(CONN);

    #Get components
    $slot=3;
    open(COM,$components) || die "Cannot open COM at $components.\n";
    while(<COM>){
      if($_=~/^$TF\t(\d+)\t/){
	$data{$slot}=$1;
	last;
      }
    }
    close(COM);
      
    #Calculate corrected connectivity
    $slot=4;
    $total=$data{1} + $data{3};
    $total--; #Only consider extra components, 1 is not relevant and onlt eliminates cnnectivities of 1.0.  
    if(!($total)){
      $corrected=0;
    }else{
      $corrected=$data{2} / $total;
    }
    #Round decimals
    $rounded = sprintf "%.2f", $corrected;
    #Save variable
    $data{$slot}=$rounded;

    #Get total catalysis reactions
    $slot=5;
    open(CAT,$reactions) || die "Cannot open CAT at $reactions\n";
    while(<CAT>){
      if($_=~/^$TF\t(\d+)/){
	$data{$slot}=$1;
	last;
      }
    }
    close(CAT);
    
    #Get secondary reactions
    $slot=6;
    open(SRX,$srs) || die "Cannot open SRX at $srs\n";
    while(<SRX>){
      if($_=~/^$TF\t(\d+)/){
	$data{$slot}=$1;
	last;
      }
    }
    close(SRX);

    #Find feedback
    $slot=7;
    open(FBK,$feedback) || die "Cannot open FBK at $feedback.\n";
    while(<FBK>){
      if($_=~/^$TF\t/){
	$data{$slot}=1;
      }
    }
    close(FBK);

    #Print output
    print OUT"$TF1";
    for($i=1;$i<=7;$i++){
      if($data{$i}){
	print OUT"\t$data{$i}";
      }else{
	print OUT"\t0";
      }
    }
    print OUT"\n";

	
  }
}
close(IN);

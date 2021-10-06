#################################################################################
###### 26-09-16                                                            ######
###### Connectivity Calculator with Effect                                 ######
###### Calculates connectivity from output files (avoids spreadsheet)      ######
###### Arguments:                                                          ######
###### [0] - calculations output directory                                 ######
###### [1] - Output file path                                              ######
###### Output:                                                             ######
######   f1 > TF                                                           ######
######   f2 > total reactions activated                                    ######
######   f3 > components activated                                         ######
######   f4 > number of enzymes activated                                  ######
######   f5 > activated enzymes with connectivity                          ######
######   f6 > corrected connectivity for activation                        ######
######   f7 > total reactions repressed                                    ######
######   f8 > components repressed                                         ######
######   f9 > number of enzymes repressed                                  ######
######   f10 > repressed enzymes with connectivity                         ######
######   f11 > corrected connectivity for repression                       ######
###### Notes:                                                              ######
###### History:                                                            ######
#################################################################################

$time=localtime();

#Get ARGS
$dir_path=$ARGV[0];
$outfile=$ARGV[1];
chomp($outfile);

#Open outfile
open(OUT,">$outfile") || die "Cannot open OUT at $outfile.\n";
print OUT"# From $0 on $time\n# Input = $dir_path\n# f1 > TF\n# f2 > Total activated reactions\n# f3 > Total components activated\n# f4 > Total activated enzymes\n# f5 > Total activated enzymes with connectivity\n# f6 > Corrected connectivity for activation\n# f7 > Total repressed reactions\n# f8 > Total components repressed\n# f9 > Total repressed enzymes\n# f10 > Total repressed enzymes with connectivity\n# f11 > Corrected connectivity for repression\n";

#Get files
$ac_components=$dir_path . "activation_components.txt";
$rs_components=$dir_path . "inhibition_components.txt";
$enzyme_count=$dir_path . "enzyme_counts.txt";
$ac_conn=$dir_path . "connectivity_activation.txt";
$rs_conn=$dir_path . "connectivity_inhibition.txt";

open(EZC,$enzyme_count) || die "Cannot open EZC at $enzyme_count.\n";
while(<EZC>){
  if($_=~/^([^\t]+)\t(\d+)\t(\d+)\t[^\t]*\t(\d+)\t(\d+)\t[^\t]*/){
    $TF=$1;
    $en_ac=$2;
    $en_rs=$4;
    $t_ac=$3;
    $t_rs=$5;

    #Reset values
    $cn_ac=0;
    $cm_ac=0;
    $cn_rs=0;
    $cm_rs=0;

    print"$TF\n";

    #Calculus for activation
    #Enzymes with connectivity
    open(CNA,$ac_conn) || die "Cannot open CNA at $ac_conn.\n";
    while(<CNA>){
      if($_=~/^$TF\t(\d+)\t/){
	$cn_ac=$1;
	last;
      }
    }
    close(CNA);

    #Components
    open(CMA,$ac_components) || die "Cannot open CMA at $ac_conmponents.\n";
    while(<CMA>){
      if($_=~/^$TF\t(\d+)/){
	$cm_ac=$1;
	last;
      }
    }
    close(CMA);

    #Avoid division by 0
    $denominator=$t_ac + ($cm_ac - 1);
    if($denominator>0){
      
      $corr_conn_ac=$cn_ac / $denominator;
      $corr_ac=sprintf "%.2f", $corr_conn_ac;
    }else{
      $corr_ac="0.00";
    }

    #Calculus for repression
    #Enzymes with connectivity
    open(CNR,$rs_conn) || die "Cannot open CNR at $rs_conn.\n";
    while(<CNR>){
      if($_=~/^$TF\t(\d+)\t/){
	$cn_rs=$1;
	last;
      }
    }
    close(CNR);

    #Components
    open(CMR,$rs_components) || die "Cannot open CMR at $rs_conmponents.\n";
    while(<CMR>){
      if($_=~/^$TF\t(\d+)\t/){
	$cm_rs=$1;
	last;
      }
    }
    close(CMR);

    #Avoid division by 0
    $denominator=$t_rs + ($cm_rs - 1);
    if($denominator>0){
      #print"den: $t_rs + ($cm_rs - 1) = $denominator\n";
      $corr_conn_rs=$cn_rs / $denominator;
      $corr_rs=sprintf "%.2f", $corr_conn_rs;
    }else{
      $corr_rs="0.00";
    }

    print OUT"$TF\t$en_ac\t$cm_ac\t$t_ac\t$cn_ac\t$corr_ac\t$en_rs\t$cm_rs\t$t_rs\t$cn_rs\t$corr_rs\n";

  }
}
close(EZC);
close(OUT);


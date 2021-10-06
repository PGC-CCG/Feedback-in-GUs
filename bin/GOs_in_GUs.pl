#################################################################################
###### 29-09-16                                                            ######
###### GO coverage per GU                                                  ######
###### Calculates the percentage of GO genes present in a GU               ######
###### Arguments:                                                          ######
###### [0]> GOs_table.txt                                                  ######
###### [1]> TF-genes RegulonDB dataset                                     ######
###### [2]> Output file                                                    ######
###### Output:                                                             ######
###### [1]> TF                                                             ######
###### [2]> GO ID & name                                                   ######
###### [3]> Genes associated to GO in Regulon                              ######
###### [4]> Percentage of GO covered                                       ######
###### History:                                                            ######
#################################################################################

$time=localtime();

#Get args
$table=$ARGV[0];
$dataset=$ARGV[1];
$outfile=$ARGV[2];
chomp($outfile);

#Create and print header of OUT file
open(OUT,">$outfile") || die "Cannot open OUT file at $outfile\n";
print OUT"# From $0 on $time\n# Input = $table / $dataset \nGU\tGO ID\tGO name\tFraction of GO covered by GU\tSize of GO\tGO-GU Overlap\tGenes in overlap\n";

undef %TFs;

#Retrieve regulons
open(TRN,$dataset) || die "Cannot open TRN at $dataset.\n";
while(<TRN>){
  #if($_=~/^([^\t]+)\t([^\t]+)\t/){
  if($_=~/^(\w+)\t([\w]+)/){
    $TF=$1;
    
    $flag=0;
    foreach $tfs (keys %TFs){
      if($tfs eq $TF){
	$flag=1;
	last;
      }
    }
    if($flag){
      next;
    }else{
      print"$TF\n";
      $TFs{$TF}++;
    }

    undef %regulon;
    undef %gos;

    #Save all regulated genes
    open(TRN1,$dataset) || die "Cannot open TRN1 at $dataset.\n";
    while(<TRN1>){
      if($_=~/^$TF\t([\w]+)/){
	$regulon{$1}++;
      }
    }
    close(TRN1);

    #Find GOs of genes
    foreach $gn (keys %regulon){

      open(TBL,$table) || die "Cannot open TBL at $table.\n";
      while(<TBL>){
	if($_=~/^([^\t]+)\t([^\t]+)\t([^\t]+)\t(.*,*$gn,[^\t]*)\t/){
	  $go=$1;
	  $go_name=$3;
	  $genes=$4;

	  #Avoid repetition of GO annotations
	  $flag=0;
	  foreach $rp_go (keys %gos){
	    if($rp_go eq $go){
	      $flag=1;
	    }
	  }
	  
	  if($flag){
	    next;
	  }else{
	    $gos{$go}++;
	  }

	  @gns=split(/,/,$genes);
	  
	  #Get percentage of GO included in regulon
	  $size=@gns;
	  @overlap=();
	  foreach $gogn (@gns){
	    foreach $rgn (keys %regulon){
	      if($rgn eq $gogn){
		push(@overlap,$rgn);
	      }
	    }
	  }
	  
	  $s_overlap=@overlap;
	  $fraction=$s_overlap/$size;
	  $percentage=sprintf "%.2f", $fraction;

	  print OUT"$TF\t$go\t$go_name\t$percentage\t$size\t$s_overlap\t";
	  foreach $ov (@overlap){
	    print OUT"$ov,";
	  }
	  print OUT"\n"; 
	}
      }
      close(TBL);

    }
  }
}
close(TRN);
close(OUT);


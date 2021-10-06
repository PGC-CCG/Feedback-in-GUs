#################################################################################
###### 26-09-16                                                            ######
###### Enzyme Count Effect                                                 ######
###### Get enzyme count considering type of regulation (+/-).              ######
###### [0] - gu_assembly_vX.pl output directory                            ######
###### [1] - Output file path                                              ######
###### Output:                                                             ######
######   f1 > TF                                                           ######
######   f2 > total activated reactions                                    ######
######   f3 > number of enzymes activated                                  ######
######   f4 > enzymes activated                                            ######
######   f5 > total activated reactions                                    ######
######   f6 > number of enzymes inhibited                                  ######
######   f7 > enzymes inhibited                                            ######
###### Notes:                                                              ######
###### - Complexes are considered as one enzyme                            ######
###### - Complexes with both effects are counted in both lists.            ######
###### - Enzymes under dual regulatios are in both lists.                  ######
###### History:                                                            ######
#################################################################################

$time=localtime();

#Get args
$dir_path=$ARGV[0];
$outfile=$ARGV[1];
chomp($outfile);

#Open outfile
open(OUT,">$outfile") || die "Cannot open OUT at $outfile.\n";
print OUT"# From $0 on $time\n# Input = $dir_path\n# f1 > TF\n# f2 > Total activated reactions\n# f3 > Total activated enzymes\n# f4 > Activated enzymes\n# f5 > Total repressed reactions\n# f6 > Total repressed enzymes\n# f7 > Repressed enzymes\n";

#Open dir 
opendir(DIR,$dir_path) || die "Cannot open DIR at $in_dir.\n";
@GU_files=readdir(DIR);
foreach $file (@GU_files){
  if(!($file=~/^\./)){

    print"\n**$file\n";

    undef %activation;
    undef %repression;
    $c_ac=0;
    $c_rs=0;
    
    #Get file names from GU folder
    $products=$dir_path . $file . "/reactants_products.txt";
    $objects=$dir_path . $file . "/objects.txt";
    $reactions=$dir_path . $file . "/reactions.txt";
    $modifications=$dir_path . $file . "/modification.txt";

    $TF=$file;

    #Find catalysis reactions
    open(RXN,$reactions) || die "Cannot open RXN at $reactions.\n";
    while(<RXN>){
      #Save reversibility and filter STATE_TRANSITION reactions
      if($_=~/^(re\d+)\t(STATE_TRANSITION|TRANSPORT)\t/i) {
	$rex=$1;
	
	#Find catalyzer of reaction
	open(MOD,$modifications) || die "Cannot open MOD at $modifications.\n";
	while(<MOD>){
	  if($_=~/^$rex\tCATALYSIS\t([^\t]+)/){
	    $enzyme=$1;
	    chomp($enzyme);
	    

	    #Discern proteins from protein complexes
	    if($enzyme=~/\w-\w/){
	      #print"ENzyme complex: $enzyme\n";
	      @monomers=split(/-/,$enzyme);
	      #Send monomers to subfunction
	      $cpx_activated=0;
	      $cpx_inhibited=0;
	      foreach $mn (@monomers){
		$effect=&reg_effect($rex,$mn);
		if($effect=~/^PHYSICAL_STIMULATION|DUAL/){
		  $cpx_activated=1;
		}
		if($effect=~/^INHIBITION|DUAL/){
		  $cpx_inhibited=1;
		}
	      }

	      #Count complexes
	      #Activation
	      if($cpx_activated){
		$activation{$enzyme}++;
		$c_ac++;
	      }
	      if($cpx_inhibited){
		$repression{$enzyme}++;
		$c_rs++;
	      }
	    }else{ #Enzyme not a complex
	      $effect=&reg_effect($rex,$enzyme);
   
	      if($effect=~/^PHYSICAL_STIMULATION|DUAL/){
		$activation{$enzyme}++;
		$c_ac++;
	      }
	      if($effect=~/^INHIBITION|DUAL/){
		$repression{$enzyme}++;
		$c_rs++;
	      }
	    }
	  }
	}
	close(MOD);
      }
    }
    close(RXN);

    #Print hashes
    @keys=keys %activation;
    $size=@keys;
    print OUT"$TF\t$c_ac\t$size\t$keys[0]";
    for($i=1;$i<$size;$i++){
      print OUT",$keys[$i]";
    }
    @keys=keys %repression;
    $size=@keys;
    print OUT"\t$c_rs\t$size\t$keys[0]";
    for($i=1;$i<$size;$i++){
      print OUT",$keys[$i]";
    }
    print OUT"\n";
    
  }
}
closedir(DIR);
close(OUT);


###############################
### Subroutine get protein regulation
### Navigates through GU_assembly files and returns regulation
### input: 
### [0] > $rex
### [1] > protein (monomer)
### output: inhibition/activation/dual 
###############################

sub reg_effect{
  
  local($rex,$enzyme,$enz_rxn,$mRNA,$trans_rxn,$effect,$activation,$repression)=($_[0],$_[1],"","","","",0,0);
  
  #Find TU of catalyzer to search transcription reaction
  open(RPD,$products) || die "Cannot open RPD at $products.\n";
  while(<RPD>){

    #Find translation reaction
    if($_=~/^([^\t]+)\t$enzyme\tproduct$/){
      $enz_rxn=$1;
      
      #Find TU (reactant of translation reaction)
      open(RPD1,$products) || die "Cannot open RPD1 at $products.\n";
      while(<RPD1>){
	if($_=~/^$enz_rxn\t([^\t]+)\treactant$/){
	  $mRNA=$1;
	  close(RPD1);
	}
      }
      #Find transcription reaction
      open(RPD2,$products) || die "Cannot open RPD2 at $products.\n";
      while(<RPD2>){
	if($_=~/^([^\t]+)\t$mRNA\tproduct$/){
	  $trans_rxn=$1;
	  close(RPD2);
	}
      }
      #Find regulatory effect (activation/repression)
      open(MOD1,$modifications) || die "Cannot open MOD1 at $modifications.\n";
      while(<MOD1>){
	if($_=~/^$trans_rxn\t([^\t]+)\t/){
	  $effect=$1;
	  if($effect=~/^PHYSICAL_STIMULATION/){
	    $activation=1;
	  }elsif($effect=~/^INHIBITION/){
	    $repression=1;
	  }elsif($effect=~/^UNKNOWN/){
	    $activation=1;
	    $repression=1;
	  }
	}
      }
      close(MOD1);

    }
  }
  close(RPD);

  if(($activation) && (!($repression))){
    return("PHYSICAL_STIMULATION");
  }elsif((!($activation)) && ($repression)){
    return("INHIBITION");
  }elsif(($activation) && ($repression)){
    return("DUAL");
  }else{
    print"REPORT: $enzyme in $rex is not regulated by $TF\n";
  }
}

#################################################################################
###### 03-11-15                                                            ######
###### Get all short paths + effect                                        ######
###### Per GU finds short-paths > X                                        ######
###### [0] - gu_assembly_vX.pl output directory                            ######
###### [1] - Output file path                                              ######
###### [2] - Path minimum length (use "2" for enhanced predictions,        ######
######       "0" for components search to include disconnected rxns)       ######
###### [3] - Effect to print: +|-|all                                      ######
###### History:                                                            ######
###### v1/05-11-15/> Rebuilt subroutines to reactant/product               ######
###### v2/09-11-15/> Output file path as:                                  ######
######               [1] > reaction path                                   ######
######               [2] > metabolite names path                           ######
######               [3] > metabolite ids path                             ######
###### v3/29-03-16/> Considered reversibility of reactions                 ######
######             > Minimum path length as input                          ######
######             > TF-effector complex formation reactions are omitted   ######
###### v4/20-09-16/> Only considers short paths that include genes under   ######
######               the same TF effect (repression/activation). Dual and  ######
######               unknown regulation are always permitted.              ######
#################################################################################

$time=localtime();

#Get args
$dir_path=$ARGV[0];
$out_file=$ARGV[1];
$limit=$ARGV[2];
$p_effect=$ARGV[3];
chomp($p_effect);

@temp_files=();

#Get temp file names (for reversible reactions and regulatory effect)
if($out_file=~/^(.+)\.txt/){
  $temp_files[0]=$1 . "_temp_act.txt";
  $temp_files[1]=$1 . "_temp_rep.txt";
}else{
  die "Cannot create temp file, outfile is not .txt\n";
}

#Open outfile
open(OUT,">$out_file") || die "Cannot open OUT at $out_file.\n";
print OUT"# From $0 on $time\n# Input = $dir_path\n# Limit of short paths = $limit\n";

#Open dir 
opendir(DIR,$dir_path) || die "Cannot open DIR at $in_dir.\n";
@GU_files=readdir(DIR);
foreach $file (@GU_files){
  if(!($file=~/^\./)){

    print"\n**$file\n";

    undef %paths;
    undef %mol_paths;
    
    #Get file names from GU folder
    $products=$dir_path . $file . "/reactants_products.txt";
    $objects=$dir_path . $file . "/objects.txt";
    $reactions=$dir_path . $file . "/reactions.txt";
    $modifications=$dir_path . $file . "/modification.txt";

    $TF=$file;

    #Create temp files
    open(TMA,">$temp_files[0]") || die "Cannot open TMA file at $temp_files[0].\n";
    open(TMR,">$temp_files[1]") || die "Cannot open TMR file at $temp_files[1].\n";

    undef %r_activation;
    undef %r_repression;

    #Create temp reactants_products file with duplicated entries for REV reactions and common regulation (repression/activation).
    open(RXN,$reactions) || die "Cannot open RXN at $reactions.\n";
    while(<RXN>){
      #Save reversibility and filter STATE_TRANSITION reactions
      if($_=~/^(re\d+)\t(STATE_TRANSITION|TRANSPORT)\tRVB\t/i) {
	$rex=$1;
	$rvb=1;
      }elsif($_=~/^(re\d+)\t(STATE_TRANSITION|TRANSPORT)\tL2R\t/){
	$rex=$1;
	$rvb=0;
      }elsif($_=~/^(re\d+)\tSUPER_REACTION\t/){
	$rex=$1;
	#Print super-reactions in both outfiles.
	open(RPDS,$products) || die "Cannot open RPDS at $products.\n";
	while(<RPDS>){
	  if($_=~/^$rex\t/){
	    if($p_effect eq "+"){
	      print TMA"$_";
	    }elsif($p_effect eq "-"){
	      print TMR"$_";
	    }elsif($p_effect eq "all"){
	      print TMA"$_";
	      print TMR"$_";
	    }else{
	      die "No effect was specified. Please specify:\n+ for activation\n- for inhibition\nall for both activation and inhibition\n";
	    }
	  }
	}
	close(RPDS);
      }
	
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
	      }elsif($effect=~/^INHIBITION|DUAL/){
		$cpx_inhibited=1;
	      }
	    }

	    #Print reactions regulated by complexes
	    #Activation output
	    if(($cpx_activated) && (!($r_activation{$rex}))){
	      #Print RPD 
	      if($rvb){
		open(RPD3,$products) || die "Cannot open RPD3 at $products.\n";
		while(<RPD3>){
		  if($_=~/^$rex\t([^\t]+)\treactant/){
		    if($p_effect=~/all|\+/){
		      print TMA"$_$rex\_rv\t$1\tproduct\n";
		      $r_activation{$rex}++;
		    }
		  }elsif($_=~/^$rex\t([^\t]+)\tproduct/){
		    if($p_effect=~/all|\+/){
		      print TMA"$rex\_rv\t$1\treactant\n$_";
		    }
		  }
		}
		close(RPD3);
	      }else{
		open(RPD3,$products) || die "Cannot open RPD3 at $products.\n";
		while(<RPD3>){
		  if($_=~/^$rex\t/){
		    if($p_effect=~/all|\+/){
		      print TMA"$_";
		      $r_activation{$rex}++;
		    }
		  }
		}
		close(RPD3);
	      }
	    }
	    #Inhibition output. Duals will be printed in both documents.
	    if(($cpx_inhibited) && (!($r_repression{$rex}))){
	      #Print RPD 
	      if($rvb){
		open(RPD3,$products) || die "Cannot open RPD3 at $products.\n";
		while(<RPD3>){
		  if($_=~/^$rex\t([^\t]+)\treactant/){
		    if($p_effect=~/all|\-/){
		      print TMR"$_$rex\_rv\t$1\tproduct\n";
		      $r_repression{$rex}++;
		    }
		  }elsif($_=~/^$rex\t([^\t]+)\tproduct/){
		    if($p_effect=~/all|\-/){
		      print TMR"$rex\_rv\t$1\treactant\n$_";
		    }
		  }
		}
		  close(RPD3);
	      }else{
		open(RPD3,$products) || die "Cannot open RPD3 at $products.\n";
		while(<RPD3>){
		  if($_=~/^$rex\t/){
		    if($p_effect=~/all|\-/){
		      print TMR"$_";
		      $r_repression{$rex}++;
		    }
		  }
		}
		close(RPD3);
	      }
	    }	  
	  }else{
	    #Homomeric complexes
	    $effect=&reg_effect($rex,$enzyme);
	    if(($effect=~/^PHYSICAL_STIMULATION|DUAL/) && (!($r_activation{$rex}))){
	      #Print RPD 
	      if($rvb){
		open(RPD3,$products) || die "Cannot open RPD3 at $products.\n";
		while(<RPD3>){
		  if($_=~/^$rex\t([^\t]+)\treactant/){
		    if($p_effect=~/all|\+/){
		      print TMA"$_$rex\_rv\t$1\tproduct\n";
		      $r_activation{$rex}++;
		    }
		  }elsif($_=~/^$rex\t([^\t]+)\tproduct/){
		    if($p_effect=~/all|\+/){
		      print TMA"$rex\_rv\t$1\treactant\n$_";
		    }
		  }
		}
		close(RPD3);
	      }else{
		open(RPD3,$products) || die "Cannot open RPD3 at $products.\n";
		while(<RPD3>){
		  if($_=~/^$rex\t/){
		    if($p_effect=~/all|\+/){
		      print TMA"$_";
		      $r_activation{$rex}++;
		    }
		  }
		}
		close(RPD3);
	      }
	    }
	    if(($effect=~/^INHIBITION|DUAL/) && (!($r_repression{$rex}))){
	      #Print RPD 
	      if($rvb){
		open(RPD3,$products) || die "Cannot open RPD3 at $products.\n";
		while(<RPD3>){
		  if($_=~/^$rex\t([^\t]+)\treactant/){
		    if($p_effect=~/all|\-/){
		      print TMR"$_$rex\_rv\t$1\tproduct\n";
		      $r_repression{$rex}++;
		    }
		  }elsif($_=~/^$rex\t([^\t]+)\tproduct/){
		    if($p_effect=~/all|\-/){
		      print TMR"$rex\_rv\t$1\treactant\n$_";
		    }
		  }
		}
		  close(RPD3);
	      }else{
		open(RPD3,$products) || die "Cannot open RPD3 at $products.\n";
		while(<RPD3>){
		  if($_=~/^$rex\t/){
		    if($p_effect=~/all|\-/){
		      print TMR"$_";
		      $r_repression{$rex}++;
		    }
		  }
		}
		close(RPD3);
	      }
	    }
	  }
	}
      }
      close(MOD);
    }
    close(RXN);
    close(TMA);
    close(TMR);
  
    #Repeat search for each TMP file (inhibition and activation).
    for($ix=0;$ix<2;$ix++){   #For instead of foreach to know which file is the loop in
      $temp=$temp_files[$ix];
      
      #Skip file if empty
      if($ix){
	@kys=keys %r_repression; 
	if(!(@kys)){
	  print"No repression in $TF.\n";
	  next;
	}
      }else{
	@kys=keys %r_activation; 
	if(!(@kys)){
	  print"No activation in $TF.\n";
	  next;
	}
      }

      #Find path from each product in the GU
      #Open RPD
      open(TMP,$temp) || die "Cannot open TMP at $temp.\n";
      while(<TMP>){

	if($_=~/^(re\d+)\t([^\t]+)\treactant/i){       #Find reactants
	  $rxn=$1;
	  $mol=$2;
	  $q_mol=quotemeta($mol);
	  #Eliminate objects that are not small molecules
	  open(OBJ,$objects) || die "Cannot open OBJ at $objects.\n";
	  while(<OBJ>){
	    if($_=~/^$q_mol\tSIMPLE_MOLECULE/){             #use only state_transition reactions
	      close(OBJ);
	      #print"$rxn,$mol,reactant,$rxn,$mol\n";
	      &search_path($rxn,$mol,"reactant",$rxn,$mol);    #send to sub
	    }
	  }
	  close(OBJ);
	}
      }
      close(TMP);
      
      ##Eliminate overlapping paths;keep the longest
      #In reaction hash
      foreach $k (keys %paths){
	foreach $y (keys %paths){
	  if(($k=~/^$y\/\/.+/) || ($k=~/^.+\/\/$y/)){
	    delete $paths{$y};
	  }
	}
      }
      #In molecules hash
      foreach $m (keys %mol_paths){
	$size=(@{$mol_paths{$m}});
	for($k=0;$k<$size;$k++){
	  for($y=0;$y<$size;$y++){
	    $q_y=quotemeta($mol_paths{$m}[$y]);
	    if(($mol_paths{$m}[$k]=~/^$q_y\/\/.+/) || ($mol_paths{$m}[$k]=~/^.+\/\/$q_y/)){
	      splice(@{$mol_paths{$m}},$y,1);
	    }
	  }
	}
      }
      
      ##Print remaining short paths
      foreach $k (keys %paths){
	print OUT"$file\t$k\t";
	##Print molecule paths
	foreach $m (@{$mol_paths{$k}}){
	  print OUT"$m;;";
	}
	print OUT"\t";
	foreach $m (@{$mol_paths{$k}}){
	  ##Print mol IDS
	  @ids=split(/\/\//,$m);
	  foreach $i (@ids){
	    if($i){
	      $q_i=quotemeta($i);
	      open(OBJ,$objects) || die "Cannot open OBJ at $objects.\n";
	      while(<OBJ>){	
		if($_=~/^$q_i\tSIMPLE_MOLECULE\t([^\t]+)/){
		  $id=$1;
		  chomp($id);
		  print OUT"//$id";
		  last;
		}
	      }
	      close(OBJ);
	    }
	  }
	  print OUT";;";
	}
	print OUT"\n";
      }
      ##Print TFs with no short paths for easier data handling of all GUs
      @keys=keys %paths;
      $size=@keys;
      if(!($size)){
	print OUT"$file\t\t\t\n";
      }
      
    }
  }
}
closedir(DIR);
close(OUT);

system("rm $temp_files[0] $temp_files[1]");

###############################
### Subroutine search_reactant
### Finds subsequent reaction by looking for the molecule as reactant
### input: 
### [0] > reaction(only relevant for type=reactant)
### [1] > molecule to search as product (only relevant for type=product)
### [2] > type
### [3] > reaction path
### [4] > molecule path
### output: none, saves path in %paths
###############################

sub search_path{

  local($rx,$mol,$type,$path,$m_path,$flag,$prd,$cycle)=($_[0],$_[1],$_[2],$_[3],$_[4],0,"",0);
  local(*RP1);
  local(*RP2);

  #print"==$rx/$mol/$type/$path\n";

  if($type eq "reactant"){
    #Find product from same reaction
    open(RP1,$temp) || die "Cannot open RP1 at $temp in &search_path:\n $rx, $mol, $path in reactant.\n";
    while(<RP1>){
      if($_=~/^$rx\t([^\t]+)\tproduct$/i){
	$prd=$1;
	$t_path=$m_path . "//" . $prd;       #add product to molecule path
	#print"+R:$rx,$prd,product,$path,$t_path\n";
	&search_path($rx,$prd,"product",$path,$t_path);         #call sub for product
      }
    }
    close(RP1);
  }

  if($type eq "product"){
    #Eliminate special characters
    $q_mol=quotemeta($mol);
    #Find next reaction
    open(RP2,$temp) || die "Cannot open RP2 at $temp in &search_path:\n n $rx, $mol, $path in product\n.";
    while(<RP2>){
      if($_=~/^(re\d+)\t$q_mol\treactant/){
	$rxn=$1;
	#Ignore reactants of the same reaction and same reaction in reversible order.
	if(($rx eq $rxn) || ($rxn=~/^$rx\_rv/) || ($rx=~/^$rxn\_rv/)){
	  next;
	}else{
	  #Ignore cycles:check if reaction, or its reversible version, was already on the path
	  @s_path=split(/\/\//,$path);
	  $cycle=0;              # flag for cycle in path
	  foreach $cy (@s_path){
	    if($cy eq $rxn){
 	      $cycle=1;
	      last;
	    }elsif(($cy=~/^$rxn\_rv/) || ($rxn=~/^$cy\_rv/)){
	      $cycle=1;
	      last;
	    }
	  }
	  if(!($cycle)){
	    $tpath=$path . "//" . $rxn;
	    $flag=1;
	    #print"=P:$rxn,$mol,reactant,$tpath,$m_path\n";
	    &search_path($rxn,$mol,"reactant",$tpath,$m_path);
	  }
	}
      }
    }
    close(RP2);

    if(!($flag)){
      @s_path=split(/\/\//,$path);
      $size=@s_path;
      #Ignore single reactions of TF-effector in case limit < 1.
      if(!($limit)){                 #If limit=0 and all rxns will be printed
	$flag1=0;
	@m_path=split(/\/\//,$m_path);    #split $m_path
	$m_size=@m_path;                  
	if($m_size == 2){                 #identify single reactions
	  $last=pop(@m_path);             #get product of single reaction
	  if($last=~/^$TF-/){         #if single reaction is a TF-effector complex formation
	    $flag1=1;                #turn on flag to omit
	  }
	}
	if(!($flag1)){                 #only print if flag is off...          
	  $paths{$path}++;
	  push(@{$mol_paths{$path}},$m_path);
	  return();
	}	
      }elsif($size>$limit){          #if limit>0 TF-effector rxns will only be printed when part of a path
	$paths{$path}++;
	push(@{$mol_paths{$path}},$m_path);
	return();
      }	
    }
  }
}



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
    return("OMIT");
  }else{
    print"REPORT: $enzyme in $rex is not regulated by $TF\n";
  }
}

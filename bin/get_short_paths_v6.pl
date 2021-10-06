#################################################################################
###### 03-11-15                                                            ######
###### Get all short paths                                                 ######
###### Per GU finds short-paths > X                                        ######
###### [0] - gu_assembly_vX.pl output directory                            ######
###### [1] - Output file path                                              ######
###### [2] - Path minimum length (use "1" for connectivity,                ######
######       "0" for components search to include disconnected rxns)       ######
###### History:                                                            ######
###### v1/05-11-15/> Rebuilt subroutines to reactant/product               ######
###### v2/09-11-15/> Output file path as:                                  ######
######               [1] > reaction path                                   ######
######               [2] > metabolite names path                           ######
######               [3] > metabolite ids path                             ######
###### v3/29-03-16/> Considered reversibility of reactions                 ######
######             > Minimum path length as input                          ######
######             > TF-effector complex formation reactions are omitted   ######
###### v4/31-01-17/> Bug fixed in reversible reactions                     ######
#################################################################################

$time=localtime();
print"Begin: $time\n";

#Get args
$dir_path=$ARGV[0];
$out_file=$ARGV[1];
$limit=$ARGV[2];
chomp($limit);

#Get temp file name (for reversible reactions)
if($out_file=~/^(.+)\.txt/){
  $temp=$1 . "_temp.txt";
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

    $time=localtime();
    print"Started: $time\n";

    undef %paths;
    undef %mol_paths;
    
    #Get file names from GU folder
    $products=$dir_path . $file . "/reactants_products.txt";
    $objects=$dir_path . $file . "/objects.txt";
    $reactions=$dir_path . $file . "/reactions.txt";

    $TF=$file;

    #Create temp file
    open(TMP,">$temp") || die "Cannot open TMP file at TMP.\n";

    #Create temp reactants_products file with duplicated entries for REV reactions.
    open(RXN,$reactions) || die "Cannot open RXN at $reactions.\n";
    while(<RXN>){
      if(($_=~/^(re\d+)\tSTATE_TRANSITION\tRVB\t/i) || ($_=~/^(re\d+)\tTRANSPORT\tRVB\t/i)){
	$rex=$1;
	open(RPD,$products) || die "Cannot open RPD at $products.\n";
	while(<RPD>){
	  if($_=~/^$rex\t([^\t]+)\treactant$/i){
	    print TMP"$rex\t$1\treactant\n$rex\_rv\t$1\tproduct\n";
	  }
	  if($_=~/^$rex\t([^\t]+)\tproduct$/i){
	    print TMP"$rex\_rv\t$1\treactant\n$rex\t$1\tproduct\n";
	  }
	}
	close(RPD);
      }elsif(($_=~/^(re\d+)\tSTATE_TRANSITION\tL2R\t/) || ($_=~/^(re\d+)\tTRANSPORT\tL2R\t/)){
	$rex=$1;
	open(RPD,$products) || die "Cannot open RPD at $products.\n";
	while(<RPD>){
	  if($_=~/^$rex\t[^\t]+\t[^\t]+$/i){
	    print TMP"$_";
	  }
	}
	close(RPD);
      }elsif(($_=~/^(re\d+)\tSUPER_REACTION\t/) || ($_=~/^(re\d+)\tSTATE_TRANSITION$/)){   #If reaction is 2ndry reaction or complex formation. Complex formation are needed to consider in connectivity single reactions that produce the effector.
	$rex=$1;
	$resr=$rex . "_sr"; #Add extension to omit call to &search_path from 2dry reactions and complex formation reactions
	open(RPD,$products) || die "Cannot open RPD at $products.\n";
	while(<RPD>){
	  if($_=~/^$rex\t([^\t]+\t[^\t]+)$/i){
	    print TMP"$resr\t$1";
	  }
	}
	close(RPD);
      }	
    }
    close(RXN);
    close(TMP);
	

    #Find path from each product in the GU
    #Open RPD
    open(TMP,$temp) || die "Cannot open TMP at $temp.\n";
    while(<TMP>){
      if($_=~/^(re\d+)\t([^\t]+)\treactant/i){       #Find reactants
	$rxn=$1;
	$mol=$2;
	if($mol=~/^phosphate$/){
	  print"Phosphate omitted in rxn $rxn.\n";
	  next;
	}
	$q_mol=quotemeta($mol);
	#Eliminate objects that are not small molecules
	open(OBJ,$objects) || die "Cannot open OBJ at $objects.\n";
	while(<OBJ>){
	  if($_=~/^$q_mol\tSIMPLE_MOLECULE/){             #use only state_transition reactions
	    close(OBJ);
	    #print"$rxn,$mol,reactant,$rxn,$mol\n";
	    #if($rxn=~/^(re\d+)\_rv/){
	    #  $n_rxn=$1;
	    #}
	    print"Sent $rxn\n";
	    &search_path($rxn,$mol,"reactant",$rxn,$mol);    #send to sub
	  }
	}
	close(OBJ);
      }
    }
    close(TMP);
   
    print"Path Overlapping\n";


    ##Eliminate overlapping paths;keep the longest
    #In reaction hash
    @all_paths=keys %paths;
    $size_paths=@all_paths;
    print"$size_paths to calculate.\n";

    undef %real_paths;

    foreach $k (keys %paths){
      $ov_flag=1;
      foreach $y (keys %paths){
	if(($y=~/^$k\/\/.+/) || ($y=~/^.+\/\/$k/)){
	  $ov_flag=0;
	  last;
	}
      }
      if($ov_flag){
	$real_paths{$k}++;
      }
    }

    print"Molecule Overlapping\n";

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



    print"Checking Out\n";


    ##Print remaining short paths
    foreach $k (keys %real_paths){
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
closedir(DIR);
close(OUT);

system("rm $temp");

$time=localtime();
print"End: $time\n";

###############################
### subroutine search_reactant
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

  local($rx,$mol,$type,$path,$m_path,$flag,$prd,$cycle,$filehl)=($_[0],$_[1],$_[2],$_[3],$_[4],0,"",0,"FHL");
  local(*RP1);
  local(*RP2);


  #print"==$file - $rx/$mol/$type/$path\n";

  if($type eq "reactant"){
    #Find product from same reaction
    open(RP1,$temp) || die "Cannot open RP1 at $temp in &search_path:\n $rx, $mol, $path in reactant.\n";
    while(<RP1>){
      if($_=~/^$rx\t([^\t]+)\tproduct$/i){
	$prd=$1;
	if($prd=~/^phosphate$/){
	  #print"Phosphate omitted in rxn $rx.\n";
	  next;
	}
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
      if($_=~/^(re\d+\w*)\t$q_mol\treactant/){
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
	    #if($rxn=~/^(re\d+)\_rv/){
	    #  $rxn=$1;
	    #}
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





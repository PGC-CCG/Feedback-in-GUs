#################################################################################
###### 05-09-15                                                            ######
###### GU Assembly                                                         ######
###### Extracts data from Ecocyc and generates 5 relational tables.        ######
###### Argument usage:                                                     ######
###### [0] - Mode: must be "d" or "s".                                     ######
######    d--Directory mode. Build GU's for every regulon file in directory. ####
######    s--Single mode. Build GU for a single regulon file.              ######
###### [1] - TFs file (TFs for which GU will be built) || TF name if mode=s  ####
###### [2] - Output directory for relational tables                        ######
###### [3] - Library dir path. Must include at least: TF_active, TF_inactive, ###
######       gene_links, tu_TRN, tu_links and gu_compound_dictionary files ######
###### History:                                                            ######
######   v0 > Built on Ecocyc v19.0                                        ######
######   v1 > 06/10/15 > Added ID and type of object to OBJ file           ######
######      > Added reaction directions and ID on RXN file                 ######
######   v2 > 30/10/15 > Added rDB type to &IDs to include effectors in    ######
######        dictionary syntaxis                                          ######
######      > Repeated reactions are also identified by comparing RPDs     ######
######   v3 > 09/11/15 Library files format as test_regulondb_utils.pl     ######
######      > Prefer protein name similar to gene name                     ######
######      > Numbered Warnings                                            ######
######   v4 > Avoid repeated entries in OBJ from &ecocyc and &effector     ######
######      > Consider effectors that are proteins                         ######
######      > Share repeated complexes hash in &ecocyc and &effector       ######
######      > Considered reactions without direction as reversible         ######
###### Notes:                                                              ######
###### - Requires Perlcyc & Pathway tools installed                        ######
###### - Pathey Tools should be running in -api mode                       ######
#################################################################################

use perlcyc;

###Open pathway tools, create new e.coli session.
#system("ptools -api -lisp");
$cyc = perlcyc -> new("ECOLI");

###Get args
$mode=$ARGV[0];
$tf_names=$ARGV[1];
$dir_out=$ARGV[2];
$lib=$ARGV[3];
chomp($lib);

### Get library files
$tf_active=$lib . "TF_active.txt";
$tf_inactive=$lib . "TF_inactive.txt";
$gene_links=$lib . "gene_links.txt";
$dictionary=$lib . "gu_compound_dictionary.txt";
$tu_trn=$lib . "tu_TRN.txt";
$tu_links=$lib . "tu_links.txt";

##Check correct mode usage.
if((($mode eq "d") || ($mode eq "s")) && ($dir_out) && ($tf_names) && ($lib)){
  ##Directory mode.
  if($mode eq "d"){
    open(NMS, $tf_names) || die "Cannot open NMS at $tf_names.\n";
    while(<NMS>){
      if($_=~/^(\w+[^\t]*)$/){
	$name=$1;
	chomp($name);
	print"\n****$name****\n";

	@startup_v=&startup($name, $dir_out);
	$outdir=$startup_v[0];
	if(!($outdir)){
	  next;
	}
	$TF=$startup_v[1];
	&effector($TF,$outdir);
      }
    }
  }
  ##Single mode.
  if($mode eq "s"){
    print"\n****$tf_names****\n";
    @startup_v=&startup($tf_names,$dir_out);
    $outdir=$startup_v[0];
    if(!($outdir)){
      die;
    }
    $TF=$startup_v[1];
    &effector($TF,$outdir);
  }
}else{
  die "Argument usage:\n[1] - Mode: must be \"d\" or \"s\"\nd--Directory mode. Build GU's for every regulon file in directory.\ns--Single mode. Build GU for a single regulon file.\n[2] - Read directory/file.\n[3] - Working directory.\n";
}


###########################
### startup Function
###
### - Arguments: TF name & output directory.
### - Finds TF name in TRN file.
### - Creates out dir.
###########################

sub startup {

   local($tf,$dir_out)=($_[0],$_[1],$_[2]);
   
   open(TRN,$tu_trn) || die "Cannot open TFA at $tu_trn.\n";             
   $qtf=quotemeta($tf);
   while(<TRN>){
     if($_=~/^$qtf\t/i){
       ###Name output directory.
       $d_out=$dir_out . $tf . "/";
   
       ###Create output directory.
       system("mkdir $d_out");

       close(TRN);
       return($d_out,$tf);
     }
   }
   close(TRN);
   print"Warning 5! $tf is not in the database, impossible to generate GU.\nPlease double-check spelling.\n";
   return(0);
 }



###########################
### Effector Function
###
### - Looks for active and inactive conformations of the TF.
### - Prints:
###     +TF-effector complexes (active & inactive)
###     +Transcription of TUs
### Arguments: 
###        [0]> TF name
###        [1]> Output dir
### Returns:
###        Vector with gene names
########################### 

sub effector {

  local($TF,$d_out)=($_[0],$_[1]);

  ###Create output files.
  $f_objects=$d_out . "ob_temp.txt";
  $f_reactions=$d_out . "reactions.txt";
  $f_complexes=$d_out . "complexes.txt";
  $f_recprod=$d_out . "reactants_products_temp.txt";
  $f_modification=$d_out . "t_modification.txt";

  ##Open output files.
  open(OBJ,">$f_objects") || die "Cannot open OBJ file for $TF.\n";
  open(RXN,">$f_reactions") || die "Cannot open RXN file for $TF.\n";
  open(CPX,">$f_complexes") || die "Cannot open CPX file for $TF.\n";
  open(RPD,">$f_recprod") || die "Cannot open RPD file for $TF.\n";
  open(MOD,">$f_modification") || die "Cannot open MOD file for $TF.\n";

  ##Reset counters.
  $c_rxn=1;
  $c_cpx=1;
  @rp_tus=();
  @effectors=();
  undef %tus;
  undef %r_prots;
  undef %r_cpx;

  ###Look for TF in TF_active.txt (active conformations).
  open(ACT,$tf_active) || die "Cannot open ACT file at $tf_active.\n";
  while(<ACT>){
    if($_=~/^$TF\t([^\t]+)\t([^\t]*)\t([^\t]+)\t([^\t]+)/i){
      #print"$_";
      $cpx=$1;
      $eff=$2;
      $tu=$4;
      $func=$3;
      if(!($eff)){
	if($cpx=~/^$TF-(.+)$/){
	  $eff=$1;
	}
      }
      if($eff){
	if($eff=~/^\w{3,6}$/){
	  $b_prot=&is_protein($eff);
	  if($b_prot){
	    print"Effector $eff labeled as protein.\n";
	    print OBJ"$eff\tPROTEIN\n";
	    push(@{$r_prots{$eff}},"EFF");
	  }else{
	    $eff=&IDS($eff,"rDB");
	    print OBJ"$eff\tSIMPLE_MOLECULE\n";
	  }
	}else{
	  $eff=&IDS($eff,"rDB");
	  print OBJ"$eff\tSIMPLE_MOLECULE\n";
	}
	$tf_cpx=$TF . "-" . $eff;
	$rpt=0;
	foreach $rp_cpx (keys %r_cpx){
	  if($tf_cpx eq $rp_cpx){
	    $rpt=1;
	    last;
	  }
	}
	if(!($rpt)){
	  push(@effectors,$eff);
	  ##Print Complex
	  print OBJ"$TF\tPROTEIN\n$tf_cpx\tCOMPLEX\n";
	  print RXN"re$c_rxn\tSTATE_TRANSITION\n";
	  print RPD"re$c_rxn\t$TF\treactant\nre$c_rxn\t$eff\treactant\nre$c_rxn\t$tf_cpx\tproduct\n";
	  print CPX"csa$c_cpx\t$tf_cpx\t$TF\ncsa$c_cpx\t$tf_cpx\t$eff\n";
	  $c_rxn++;
	  $c_cpx++;
	  $r_cpx{$tf_cpx}++;
	  push(@{$r_prots{$TF}},"EFF");
	}
      }else{
	##Print TF as object.
	print OBJ"$TF\tPROTEIN\n";
	push(@{$r_prots{$TF}},"EFF");
	$tf_cpx=$TF;
      }

      $rpt=0;
      foreach $rp_tu (@rp_tus){
	if($tu eq $rp_tu){
	  $rpt=1;
	  last;
	}
      }
      if(!($rpt)){
	##Print TF effect over TU
	if(($func eq "+") || ($func eq "+?")){
	  print MOD"re$c_rxn\tPHYSICAL_STIMULATION\t$tf_cpx\n";
	}elsif(($func eq "-") || ($func eq "-?")){
	  print MOD"re$c_rxn\tINHIBITION\t$tf_cpx\n";
	}elsif($func eq "+-"){
	  print MOD"re$c_rxn\tINHIBITION\t$tf_cpx\n";
	  print MOD"re$c_rxn\tPHYSICAL_STIMULATION\t$tf_cpx\n";
	}elsif($func eq "?"){
	  print MOD"re$c_rxn\tUNKNOWN\t$tf_cpx\n";
	}elsif(!($func)){
	  print MOD"re$c_rxn\tUNKNOWN\t$tf_cpx\n";
	}
	print OBJ"$tu\tGENE\n$tu\_mRNA\tRNA\n";
	print RXN"re$c_rxn\tTRANSCRIPTION\n";
	print RPD"re$c_rxn\t$tu\treactant\nre$c_rxn\t$tu\_mRNA\tproduct\n";
	push(@rp_tus,$tu);
	$tus{$tu}=$c_rxn;
	$c_rxn++;
      }else{
	$rx=$tus{$tu};
	if(($func eq "+") || ($func eq "+?")){
	  print MOD"re$rx\tPHYSICAL_STIMULATION\t$tf_cpx\n";
	}elsif(($func eq "-") || ($func eq "-?")){
	  print MOD"re$rx\tINHIBITION\t$tf_cpx\n";
	}elsif($func eq "+-"){
	  print MOD"re$rx\tINHIBITION\t$tf_cpx\n";
	  print MOD"re$rx\tPHYSICAL_STIMULATION\t$tf_cpx\n";
	}elsif($func eq "?"){
	  print MOD"re$rx\tUNKNOWN\t$tf_cpx\n";
	}elsif(!($func)){
	  print MOD"re$rx\tUNKNOWN\t$tf_cpx\n";
	}
      }
    }
  }
  close(ACT);
  
  ##Print inactive complexes
  open(ITV,$tf_inactive) || die "Cannot open ITV file at $tf_inactive.\n";
  while(<ITV>){
    if($_=~/^$TF\t([^\t]+)/i){
      $i_eff=$1;
      chomp($i_eff);
      $ef="";
      @inactive=split(/\|/,$i_eff);
      foreach $in_cpx (@inactive){
	#Ignore TF as complex
	if($in_cpx=~/^$TF$/i){
 	  next;
	#Get effector
	}elsif($in_cpx=~/^$TF-(.+)$/){
	  $eff=$1;
	  if($eff=~/^\w{3,6}$/){
	    
	    $b_prot=&is_protein($eff);
	    if($b_prot){
	      print"Effector $eff labeled as protein.\n";
	      print OBJ"$eff\tPROTEIN\n";
	      push(@{$r_prots{$eff}},"EFF");
	    }else{
	      $eff=&IDS($eff,"rDB");
	      print OBJ"$eff\tSIMPLE_MOLECULE\n";
	    }
	  }else{
	    $eff=&IDS($eff,"rDB");
	    print OBJ"$eff\tSIMPLE_MOLECULE\n";
	  }
	  $in_cpx=$TF . "-" . $eff;
	}else{
	  print"WARNING 1! Inactive complex $in_cpx does not follow TF-effector format and was omitted.\n";
	  next;
	}
	#Omit repeated complexes
	$rpt=0;
	foreach $rp_cpx (keys %r_cpx){
	  if($in_cpx eq $rp_cpx){
	    $rpt=1;
	    last;
	  }
	}
	if(!($rpt)){
	  push(@effectors,$eff);
	  print OBJ"$TF\tPROTEIN\n$in_cpx\tCOMPLEX\n";
	  print RXN"re$c_rxn\tSTATE_TRANSITION\n";
	  print RPD"re$c_rxn\t$TF\treactant\nre$c_rxn\t$eff\treactant\nre$c_rxn\t$in_cpx\tproduct\n";
	  print CPX"csa$c_cpx\t$in_cpx\t$TF\ncsa$c_cpx\t$in_cpx\t$eff\n";
	  $c_rxn++;
	  $c_cpx++;
	  $r_cpx{$in_cpx}++;
	  push(@{$r_prots{$TF}},"EFF");
	}
      }
    }
  }
  close(ITV);

  ##Get regulatory interactions without effector.
  open(TRN,$tu_trn) || die "Cannot open TRN at $tu_trn.\n";
  while(<TRN>){
    if($_=~/^$TF\t([^\t]+)\[[^\t]+\t([^\t]+)/i){
      $tu=$1;
      $func=$2;
      $rpt=0;
      foreach $rt (@rp_tus){
	if($tu eq $rt){
	  $rpt=1;
	  last;
	}
      }
      $tf_cpx=$TF;
      if(!($rpt)){
	##Print TF effect over TU
	if($func eq "+"){
	  print MOD"re$c_rxn\tPHYSICAL_STIMULATION\t$tf_cpx\n";
	}elsif($func eq "-"){
	  print MOD"re$c_rxn\tINHIBITION\t$tf_cpx\n";
	}elsif($func eq "+-"){
	  print MOD"re$c_rxn\tINHIBITION\t$tf_cpx\n";
	  print MOD"re$c_rxn\tPHYSICAL_STIMULATION\t$tf_cpx\n";
	}elsif($func eq "?"){
	  print MOD"re$c_rxn\tUNKNOWN\t$tf_cpx\n";
	}elsif(!($func)){
	  print MOD"re$c_rxn\tUNKNOWN\t$tf_cpx\n";
	}
	print OBJ"$tu\tGENE\n$tu\_mRNA\tRNA\n";
	print RXN"re$c_rxn\tTRANSCRIPTION\n";
	print RPD"re$c_rxn\t$tu\treactant\nre$c_rxn\t$tu\_mRNA\tproduct\n";
	push(@rp_tus,$tu);
	$tus{$tu}=$c_rxn;
	$c_rxn++;
      }
    }
  }
  close(TRN);
  close(OBJ);
  close(RXN);
  close(CPX);
  close(RPD);
  close(MOD);


  $done=&ecocyc(\%tus,$d_out,\%r_cpx,$c_rxn,$c_cpx,\%r_prots);
  if($done){
    return($done);
  }else{
    die "Error in Ecocyc function.\n $TF did not return positive value.\n";
  }

}

###########################
### Ecocyc Function
###
### - Looks for gene products, complexes and the reactions they catalyze.
### - Prints:
###     +Gene Products/Complexes
###     +Reactions
### Arguments: 
###        [0]> Gene vector ref 
###        [1]> Output dir
### Returns:
###        Termination value (1).
########################### 

sub ecocyc {

  local(%tus)=(%{$_[0]});
  local($d_out)=($_[1]);
  local(%r_cpx)=(%{$_[2]});
  local($c_rxn)=($_[3]);
  local($c_cpx)=($_[4]);
  local(%r_prots)=(%{$_[5]});

  ###Create temporal output file.
  $temp=$d_out . "data_temp.txt";
  open(TMP,">$temp") || die "Cannot open TMP file at $temp\n";

  ##Get genes
  foreach $tu (keys (%tus)){
    $qtu=quotemeta($tu);
    open(TLK,$tu_links) || die "Cannot open TLK at $tu_links.\n";
    $tu_present=0;    #Flag to recognize TUs absent from tu_links.txt
    while(<TLK>){
      if($_=~/^[^\t]+\t$qtu\t[^\t]+\t([^\t]+)\t/){
	$tu_present=1;
	$gns=$1;	
	@genes=split(/,/,$gns);
	foreach $gn (@genes){
	  $g=&IDS($gn,"g");
	  @gprods=();
	  @gprods=$cyc->all_products_of_gene($g);
	  foreach $en (@gprods){
	    #Ignore if product is a complex (will appear below from a monomer)
	    #Needed for identifying the regulated monomer
	    @mons=();
	    @mons=$cyc->monomers_of_protein($en);
	    $size=@mons;
	    #Separate homomeric from heteromultimeric complexes. Ignore hetero-.
	    if($size>1){
	      next;
	    }
	    #Print enzyme; needed for translation reaction
	    print TMP"$tu\t$en\n";
	    #Check if product is part of a complex
	    @cpxs=();
	    @cpxs=$cyc->get_slot_values($en,"COMPONENT-OF");
	    foreach $cp (@cpxs){
	      @mons=();
	      $all_mons="";
	      @mons=$cyc->monomers_of_protein($cp);
	      $size=@mons;
	      if($size>1){
		foreach $mn (@mons){
		  $all_mons=$all_mons . "$mn,";
		}
		print TMP"$tu\t$en\t$cp\t$all_mons\n";
		$ref=&enzyme($cp);
		@rxns=(@{$ref});
		foreach $rx (@rxns){
		  print TMP"$tu\t$en\t$cp\t\t$rx\n";
		}
	      }else{
		$ref=&enzyme($cp);
		@rxns=(@{$ref});
		foreach $rx (@rxns){
		  print TMP"$tu\t$en\t\t\t$rx\n";
		}
	      }
	    }
	    #Check for product reactions
	    $ref=&enzyme($en);
	    @rxns=(@{$ref});
	    foreach $rx (@rxns){
	      print TMP"$tu\t$en\t\t\t$rx\n";
	    }
	  }
	}
      }
    }
    close(TLK);
    if(!($tu_present)){
      print"WARNING!  $tu not present in tu_links for $TF \n";
    }
  }
  close(TMP);

  ###Create output files.
  $f_objects=$d_out . "ob_temp.txt";
  $f_reactions=$d_out . "reactions.txt";
  $f_complexes=$d_out . "complexes.txt";
  $f_recprod=$d_out . "reactants_products_temp.txt";
  $f_modification=$d_out . "t_modification.txt";
  $uniq_objects=$d_out . "objects.txt";
  $uniq_mods=$d_out . "modification.txt";

  ##Open output files.
  open(OBJ,">>$f_objects") || die "Cannot open OBJ file for $TF.\n";
  open(RXN,">>$f_reactions") || die "Cannot open RXN file for $TF.\n";
  open(CPX,">>$f_complexes") || die "Cannot open CPX file for $TF.\n";
  open(MOD,">>$f_modification") || die "Cannot open MOD file for $TF.\n";
  
  undef %r_rxn;

  open(TMP,$temp) || die "Cannot open TMP file at $temp\n";
  while(<TMP>){

    #Open RPD inside the loop to avoid conflicts in repeated reactions.
    open(RPD,">>$f_recprod") || die "Cannot open RPD file for $TF.\n";

    ## Print translation reactions
    if($_=~/^([^\t]+)\t([^\t]+)$/){
      $tu=$1;
      $pr_id=$2;
      chomp($pr_id);
      
      #Get protein name
      $pr_nm=&IDS($pr_id,"p");

      #Ignore repeated proteins from the same TU
      $flag=0;
      $flag1=0;
      foreach $rp_tu (@{$r_prots{$pr_nm}}){
	if($rp_tu eq $tu){
	  $flag=1;
	}
	if($rp_tu eq "EFF"){
	  $flag1=1;
	}
      }
      if($flag){
	next;
      }

      push(@{$r_prots{$pr_nm}},$tu);
      print RXN"re$c_rxn\tTRANSLATION\n";
      print RPD"re$c_rxn\t$tu\_mRNA\treactant\n";
      print RPD"re$c_rxn\t$pr_nm\tproduct\n";
      #Ignore proteins already annotated in &effectors (avoid duplicated OBJ entries)
      if(!($flag1)){
	print OBJ"$pr_nm\tPROTEIN\t$pr_id\n";
      }
      $c_rxn++;
    }

    ## Print complexes
    if($_=~/^[^\t]+\t[^\t]+\t([^\t]+)\t([^\t]+)$/){
      $cpx_id=$1;
      $monos=$2;
      chomp($monos);

      #Get complex name
      $cpx_nm=&IDS($cpx_id,"cx");
      
      #Ignore repeated complexes
      $flag=0;
      foreach $rp_cp (keys %r_cpx){
	if($rp_cp eq $cpx_nm){
	  $flag=1;
	  last;
	}
      }
      if($flag){
	next;
      }
      $r_cpx{$cpx_nm}++;

      ## Print complex formation reaction
      print RXN"re$c_rxn\tSTATE_TRANSITION\tL2R\n";
      print OBJ"$cpx_nm\tCOMPLEX\t$cpx_id\n";
      @monos=split(/,/,$monos);
      foreach $mono_id (@monos){
	$mono_nm=&IDS($mono_id,"p");

	#Omit OBJ print of proteins printed in &effectors
	$flag1=0;
	foreach $rp_pr (@{$r_prots{$pr_nm}}){
	  if($rp_pr eq "EFF"){
	    $flag1=1;
	  }
	}
	if(!($flag1)){
	  print OBJ"$mono_nm\tPROTEIN\t$mono_id\n";
	}
	print RPD"re$c_rxn\t$mono_nm\treactant\n";
	
	## Print complex id
	print CPX"csa$c_cpx\t$cpx_nm\t$mono_nm\n"
      }
      print RPD"re$c_rxn\t$cpx_nm\tproduct\n";
      $c_rxn++;
      $c_cpx++;
    }
      
    ## Print reactions
    if($_=~/^[^\t]+\t([^\t]+)\t([^\t]*)\t\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)$/){
      $mono=$1;
      $cpx=$2;
      $rxn_id=$3;
      $rxn_dir=$4;
      $reactants=$5;
      $products=$6;
      chomp($products);

      #Identify catalyzer
      if($cpx){
	$cat_id=$cpx;
	$cat_nm=&IDS($cat_id,"cx");
      }else{
	$cat_id=$mono;
	$cat_nm=&IDS($cat_id,"p");
      }

      #Ignore repeated reactions from same catalyzer 
      $flag=0;
      foreach $rp_cat (@{$r_rxn{$rxn_id}}){
	if($rp_cat eq $cat_nm){
	  $flag=1;
	  last;
	}
      }
      if($flag){
	next;
      }
      push(@{$r_rxn{$rxn_id}},$cat_nm);  

      ###Print modification of repeated reactions
      $flag=0;
      undef %rxn_index;

      ##Look for repeated reactions in RPD
      #Get substrates and products of reaction
      #Get individual substrates & products
      @substrates=split(/,/,$reactants);
      @products=split(/,/,$products);

      foreach $s (@substrates){
	if($s=~/^(.+)_Ext$/){
	  $sm=$1;
	  $s_nm=&IDS($sm,"sm");
	  $s_nm=$s_nm . "_Ext" . "//" . $s;
	}else{
	  $s_nm=&IDS($s,"sm");
	  $s_nm=$s_nm . "//" . $s;
	}
	push(@{$rxn_index{"substrates"}},$s_nm);
      }

      foreach $p (@products){
	if($p=~/^(.+)_Ext$/){
	  $sm=$1;
	  $p_nm=&IDS($sm,"sm");
	  $p_nm=$p_nm . "_Ext" . "//" . $p;
	}else{
	  $p_nm=&IDS($p,"sm");
	  $p_nm=$p_nm . "//" . $p;
	}
	push(@{$rxn_index{"products"}},$p_nm);
      }

      close(RPD);

      ##Compare to RPD reactions by sm name
      #Get RPD substrates and products by reaction
      for($inx=1;$inx<$c_rxn;$inx++){
	@rpd_subs=();
	@rpd_prods=();
	$fhl="fhl" . "$inx";
	open($fhl,$f_recprod) || die "Cannot open RPD1 at $f_recprod\n";
	while(<$fhl>){
	  if($_=~/^re$inx\t([^\t]+)\treactant/){
	    push(@rpd_subs,$1);
	  }elsif($_=~/^re$inx\t([^\t]+)\tproduct/){
	    push(@rpd_prods,$1);
	  }
	}
	close($fhl);
	
	#Compare reaction/products length
	if((@rpd_subs == @substrates) && (@rpd_prods == @products)){

	  #Compare reactants
	  $subs=@substrates;
	  $c_subs=0;
	  foreach $s (@{$rxn_index{"substrates"}}){
	    foreach $r (@rpd_subs){
	      $q_r=quotemeta($r);
	      if($s=~/^$q_r\/\/.+/i){
		$c_subs++;
	      }
	    }
	  }
	  #Compare products
	  $prods=@products;
	  $c_prods=0;
	  foreach $p (@{$rxn_index{"products"}}){
	    foreach $r (@rpd_prods){
	      $q_r=quotemeta($r);
	      if($p=~/^$q_r\/\/.+/i){
		$c_prods++;
	      }
	    }
	  }
	  
	  #Compare intersections
	  if(($c_subs == $subs) && ($c_prods == $prods)){
	    #Print modification of reaction
	    print MOD"re$inx\tCATALYSIS\t$cat_nm\n";
	    $flag=1;
	    last;
	  }
	}
      }
	 
      #If reaction has not been printed
      if(!($flag)){
	#Print reaction
	print MOD"re$c_rxn\tCATALYSIS\t$cat_nm\n";
	#Identify transport reactions
	if(($reactants=~/_Ext/) || ($products=~/_Ext/)){
	  print RXN"re$c_rxn\tTRANSPORT\t$rxn_dir\t$rxn_id\n";
	}else{
	  print RXN"re$c_rxn\tSTATE_TRANSITION\t$rxn_dir\t$rxn_id\n";
	}
	
	##Print Reactants & Products
	open(RPD,">>$f_recprod") || die "Cannot open RPD file for $TF.\n";
	foreach $s (@{$rxn_index{"substrates"}}){
	  if($s=~/^(.+)\/\/(.+)$/){
	    print OBJ"$1\tSIMPLE_MOLECULE\t$2\n";
	    print RPD"re$c_rxn\t$1\treactant\n";
	  }
	}
	foreach $p (@{$rxn_index{"products"}}){
	  if($p=~/^(.+)\/\/(.+)$/){
	    print OBJ"$1\tSIMPLE_MOLECULE\t$2\n";
	    print RPD"re$c_rxn\t$1\tproduct\n";
	  }
	}
	$c_rxn++;
	close(RPD);
      }
      

    }

  }
  close(TMP);
  close(OBJ);
  close(RXN);
  close(CPX);
  close(RPD);
  close(MOD);


  ### Sort reactants & products by length in RPD

  # Open OUT file
  $f_recprod2=$d_out . "reactants_products.txt";
  open(RPO,">$f_recprod2") || die "Cannot open RPO at $f_recprod2/\n";

  for($index=1;$index<=$c_rxn;$index++){

    # Open IN file
    open(RPD,$f_recprod) || die "Cannot open RPD file for $TF.\n";
    
    #Get objects
    @reactants=();
    @products=();
    while(<RPD>){
      if($_=~/^re$index\t([^\t]+)\treactant$/){
	push(@reactants,$1);
      }elsif($_=~/^re$index\t([^\t]+)\tproduct$/){
	push(@products,$1);
      }
    }
    close(RPD);
    
    ## Print reactants
    #Get words length
    undef %length;
    foreach $w (@reactants){
      @chars=split(//,$w);
      $size=@chars;
      $length{$w}=$size;
    }
    
    #Get longest and print in decreasing order
    @rp_word=();
    $longest="";
    $l_size=0;
    foreach $k (keys %length){
      $flag=0;
      foreach $r (@rp_word){
	if($r eq $k){
	  $flag=1;
	}
      }
      if(!($flag)){
	push(@rp_word,$k);
	if($length{$k} > $l_size){
	  $l_size=$length{$k};
	  $longest=$k;
	}
      }
    }
    
    for($i=$l_size;$i>0;$i--){
      foreach $k (keys %length){
	if($length{$k} == $i){
	  print RPO"re$index\t$k\treactant\n";
	}
      }
    }

    ## Print products
    #Get words length
    undef %length;
    foreach $w (@products){
      @chars=split(//,$w);
      $size=@chars;
      $length{$w}=$size;
    }
    
    #Get longest and print in decreasing order
    @rp_word=();
    $longest="";
    $l_size=0;
    foreach $k (keys %length){
      $flag=0;
      foreach $r (@rp_word){
	if($r eq $k){
	  $flag=1;
	}
      }
      if(!($flag)){
	push(@rp_word,$k);
	if($length{$k} > $l_size){
	  $l_size=$length{$k};
	  $longest=$k;
	}
      }
    }
    
    for($i=$l_size;$i>0;$i--){
      foreach $k (keys %length){
	if($length{$k} == $i){
	  print RPO"re$index\t$k\tproduct\n";
	}
      }
    }
  }


  #Print TF (if not printed already)
  open(OBJ,$f_objects) || die "Cannot open OBJ at $f_objects.\n";
  $flag=0;
  while(<OBJ>){
    if($_=~/^$TF\t/){
      $flag=1;
    }
  }
  close(OBJ);
  if(!($flag)){
    open(OBJ,">>$f_objects") || die "Cannot open OBJ at $f_objects.\n";
    print OBJ"$TF\tPROTEIN\n";
    close(OBJ);
  }

  #Delete duplicated objects
  system("cat $f_objects | sort | uniq > $uniq_objects");
  system("cat $f_modification | sort | uniq > $uniq_mods");
  #Remove temporal files.
  system("rm $f_objects $f_recprod $temp $f_modification");

  return(1);

}


###########################
### enzyme Function
###
### Look for enzyme reactions and annotate them
### -Input arguments:
###  [0] - enzyme name
### -Returned values: array with data_temp lines
###########################

sub enzyme{

  local($enz)=($_[0]);

  @v_out=();
  @rxns=();

  @rxns=$cyc->reactions_of_enzyme($enz);
  foreach $rx (@rxns){
    @f_phyr=();
    $all_subs="";
    $all_prods="";
    @f_phyr=$cyc->get_slot_values($rx,"REACTION-DIRECTION");
    #If no reaction direction
    if(!(@f_phyr)){
      $f_phyr[0]="REVERSIBLE";
      print"Reaction $rx without direction. Considered reversible.\n";
    }
    $f_trans=0;
    $f_trans=$cyc->transport_rxn($rx);
    if($f_trans){
      if($f_phyr[0]=~/LEFT-TO-RIGHT/){
	@subs=();
	@subs=$cyc->get_slot_values($rx,"LEFT");
	foreach $s (@subs){
	  if(!($s=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$|^\|*Pi\|*$/i)){
	    $all_subs=$all_subs . "$s\_Ext,";
	  }
	}
	@prods=();
	@prods=$cyc->get_slot_values($rx,"RIGHT");
	foreach $p (@prods){
	  if(!($p=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$|^\|*Pi\|*$/i)){
	    $all_prods=$all_prods . "$p,";
	  }
	}
	#Complete and save line
	$line="$rx\tL2R\t$all_subs\t$all_prods";
	push(@v_out,$line);
      }elsif($f_phyr[0]=~/RIGHT-TO-LEFT/){
	@subs=();
	@subs=$cyc->get_slot_values($rx,"RIGHT");
	foreach $s (@subs){
	  if(!($s=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$|^\|*Pi\|*$/i)){
	    $all_subs=$all_subs . "$s\_Ext,";
	  }
	}
	@prods=();
	@prods=$cyc->get_slot_values($rx,"LEFT");
	foreach $p (@prods){
	  if(!($p=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$|^\|*Pi\|*$/i)){
	    $all_prods=$all_prods . "$p,";
	  }
	}
	$line="$rx\tL2R\t$all_subs\t$all_prods";
	push(@v_out,$line);
      }elsif($f_phyr[0]=~/REVERSIBLE/){
	@subs=();
	@subs=$cyc->get_slot_values($rx,"LEFT");
	foreach $s (@subs){
	  if(!($s=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$|^\|*Pi\|*$/i)){
	    $all_subs=$all_subs . "$s\_Ext,";
	  }
	}
	@prods=();
	@prods=$cyc->get_slot_values($rx,"RIGHT");
	foreach $p (@prods){
	  if(!($p=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$|^\|*Pi\|*$/i)){
	    $all_prods=$all_prods . "$p,";
	  }
	}
	$line="$rx\tRVB\t$all_subs\t$all_prods";
	push(@v_out,$line);
      }
    }else{
      if($f_phyr[0]=~/LEFT-TO-RIGHT/){
	@subs=();
	@subs=$cyc->get_slot_values($rx,"LEFT");
	foreach $s (@subs){
	  if(!($s=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$/i)){
	    $all_subs=$all_subs . "$s,";
	  }
	}
	@prods=();
	@prods=$cyc->get_slot_values($rx,"RIGHT");
	foreach $p (@prods){
	  if(!($p=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$/i)){
	    $all_prods=$all_prods . "$p,";
	  }
	}
	$line="$rx\tL2R\t$all_subs\t$all_prods";
	push(@v_out,$line);
      }elsif($f_phyr[0]=~/RIGHT-TO-LEFT/){
	@subs=();
	@subs=$cyc->get_slot_values($rx,"RIGHT");
	foreach $s (@subs){
	  if(!($s=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$/i)){
	    $all_subs=$all_subs . "$s,";
	  }
	}
	@prods=();
	@prods=$cyc->get_slot_values($rx,"LEFT");
	foreach $p (@prods){
	  if(!($p=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$/i)){
	    $all_prods=$all_prods . "$p,";
	  }
	}
	$line="$rx\tL2R\t$all_subs\t$all_prods";
	push(@v_out,$line);
      }elsif($f_phyr[0]=~/REVERSIBLE/){
	@subs=();
	@subs=$cyc->get_slot_values($rx,"LEFT");
	foreach $s (@subs){
	  if(!($s=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$/i)){
	    $all_subs=$all_subs . "$s,";
	  }
	}
	@prods=();
	@prods=$cyc->get_slot_values($rx,"RIGHT");
	foreach $p (@prods){
	  if(!($p=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$/i)){
	    $all_prods=$all_prods . "$p,";
	  }
	}
	$line="$rx\tRVB\t$all_subs\t$all_prods";
	push(@v_out,$line);
      }
    }
  }

  return(\@v_out);
}



###########################
### IDS Function
###
### Turn Ecocyc id's into common names and vv.
### -Input arguments:
###  [0] - id/common name
###  [1] - type of entity (g=gene, p=protein, pn=protein to convert to gene, cx=complex, sm=simple molecule).
### -Returned values: the id or common name, opposite of input.
###########################

sub IDS {
 
  local($name,$type)=($_[0],$_[1]);

  ##Genes to ecocyc ID.
  if($type eq "g"){
    local($id);
    open(IDS, $gene_links) || die "Cannot open IDS at $gene_links.\n";
    while(<IDS>){
      if($_=~/^([^\t]+)\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t$name$/i){
	$id=$1;
	close(IDS);
	return($id);
      }
    }
    close(IDS);
    print"WARNING 2! Gene $name is not present in gene_links.\n";
    return($name);
  }

  ##Gene name from Ecocyc ID
  if($type eq "gi"){
    local($nm);
    open(IDS, $gene_links) || die "Cannot open IDS at $gene_links.\n";
    while(<IDS>){
      if($_=~/^$name\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t([^\t]+)$/i){
	$nm=$1;
	chomp($nm);
	close(IDS);
	return($nm);
      }
    }
    close(IDS);
    print"Warning 6! $name was not found at $gene_links in &IDS type \"gi\".\n";
    return($name);
  }

  ##Proteins
  if($type eq "p"){

    #RNA
    if(($name=~/RNA/) && (!($name=~/MONOMER/))){
      #Declare local variables
      local($rna_id);     
      #If class name (between pipes)
      if($name=~/^\|(.+)\|$/){
	$name=$1;
	return($name);
      }else{
	#Get synonyms
	$rna_id=$cyc->get_slot_value($name, "SYNONYMS");
	#If no synonyms
	if(!($rna_id)){
	  #Get names
	  $rna_id=$cyc->get_slot_value($name, "NAMES");
	  #If no names
	  if(!($rna_id)){
	    return("RNA");
	    #If names
	  }else{
	    $rna_id=&xml($rna_id);
	    return("RNA:$rna_id");
	  }
	  #If synonyms
	}else{
	  $rna_id=&xml($rna_id);
	  return("RNA:$rna_id");
	}
      }
    }

    ###Proteins
    #Declare local variables
    local(@gn_id,$sz_id,@p_id,$pgi,$gn_nm,$p_id);
    #Get gene to find similar name (if only 1 gene)
    @gn_id=$cyc->genes_of_protein($name);
    $sz_id=@gn_id;
    if($sz_id == 1){
      $gn_nm=&IDS($gn_id[0],"gi");
      @p_id=$cyc->get_slot_values($name,"NAMES");
      foreach $pgi (@p_id){
	if($pgi=~/^$gn_nm$/i){
	  return($pgi);
	}
      }
    }
    $p_id=$cyc->get_slot_value($name,"SYNONYMS");
    $p_id=&xml($p_id);
    #If no synonyms
    if((!($p_id=~/^\w{3,7}$/)) || ($p_id=~/^B\d+/)){
      $p_id=&IDS($name,"pn");
      return($p_id);
    }else{
      return($p_id);
    }

  }

  ##Heteromeric Complexes
  if($type eq "cx"){
    #Declare local variables
    local(@mons,$m,$c_name)=((),"","");
    @mons=$cyc->monomers_of_protein($name);
    foreach $m (@mons){
      $m_nm=&IDS($m,"p");
      $c_name=$c_name . "-" . $m_nm;
    }
    if($c_name=~/^-(.+)/){
      $c_name=$1;
      return($c_name);
    }else{
      print"Warning 7! No complex name for $name.\n";
      return"NP";
    }
  }

  ##Small Molecules from Ecocyc IDs
  if($type eq "sm"){

    #Declare local variables
    local($q_name);
    $q_name=quotemeta($name);
    open(SMS,$dictionary) || die "Cannot open SMS file at $dictionary.\n";
    while(<SMS>){
      if($_=~/^$q_name\t([^\t]+)\t/i){
	$name=$1;
	$name=&xml($name);
	return($name);
      }
    }
    close(SMS);
    if($name=~/^\|(.+)\|$/){
      return($1);
    }
    print"Warning 8! Small Molecule $name not found in dictionary.\n";
    return($name);
  }

  ##Small Molecules from Regulon DB
  if($type eq "rDB"){
    #Declare local variables
    local($q_name);
    $q_name=quotemeta($name);
    open(SMS,$dictionary) || die "Cannot open SMS file at $dictionary.\n";
    while(<SMS>){
      if($_=~/^[^\t]+\t([^\t]+)\t.*\/\/$q_name\/\//i){
	$name=$1;
	$name=&xml($name);
	return($name);
      }
    }
    close(SMS);
    print"WARNING 3! $name from RegulonDB not found in dictionary.\n";
    $name=&xml($name);
    return($name);
  }

  ##Proteins short name (gene in form XxxX)
  if($type eq "pn"){
    #Declare local variables
    local($g_pt,@g_pt,$g_nm);
    @g_pt=$cyc->genes_of_protein($name);
    open(IDS, $gene_links) || die "Cannot open IDS at $gene_links.\n";
    while(<IDS>){
      if($_=~/^$g_pt[0]\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t([^\t]+)/i){
	$g_nm=$1;
	$g_nm=ucfirst($g_nm);
	chomp($g_nm);
	return($g_nm);
      }
    }
    close(IDS);
    print"Warning 4! No gene name for $g_pt in $name.\n";
    return("NP");
  }
      

}


###########################
### XML Function
###
### Remove XML code from molecule names
### -Input arguments:
###  [0] - line
### -Returned values: line in a single char value
###########################

sub xml{

  local($line)=$_[0];
  
  chomp($line);

  while($line=~/(.*)\<.+\>(.*)/){
    $line=$1 . $2;
  }
  while($line=~/(.*)&alpha;(.*)/i){
    $line=$1 . "alpha" . $2;
  }
  while($line=~/(.*)&beta;(.*)/i){
    $line=$1 . "beta" . $2;
  }
  while($line=~/(.*)&gamma;(.*)/i){
    $line=$1 . "gamma" . $2;
  }
  while($line=~/(.*)&delta;(.*)/i){
    $line=$1 . "delta" . $2;
  }
  while($line=~/(.*)&sigma;(.*)/i){
    $line=$1 . "sigma" . $2;
  }
  while($common=~/(.*)&omega;(.*)/i){
    $common=$1 . "omega" . $2;
  }
  return($line);
  
}

###########################
### is_protein function
###
### Looks for a protein name in protein-links. Returns 0 or 1 if absent/present
### -Input arguments: protein name
### -Returned values: line in a single char value
###########################

sub is_protein {

  local($name,$q_name)=($_[0],"");

  open(PRT,$gene_links) || die "Cannot open PRT in file $gene_links at &is_protein\n";
  $q_name=quotemeta($name);
  while(<PRT>){
    if(($_=~/\t$q_name\t/i) || ($_=~/\t$q_name$/i)) {
      close(PRT);
      return(1);
    }
  }
  close(PRT);
  return(0);
}

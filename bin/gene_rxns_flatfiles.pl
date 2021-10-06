#################################################################################
###### 11-03-16                                                            ######
###### Create gene-reactions flat files                                    ######
###### Retrieves reactions of all genes, reactants and products            ######
###### Argument usage:                                                     ######
###### [0] - gene links path                                               ######
###### [1] - dictionary path                                               ######
###### [2] - outfile                                                       ######
###### Output:                                                             ######
######   [0] > gene                                                        ######
######   [1] > reaction                                                    ######
######   [2] > reactants                                                   ######
######   [3] > products                                                    ######
###### Notes:                                                              ######
###### - Requires Perlcyc & Pathway tools installed                        ######
###### - Pathway Tools should be running in -api mode                       ######
#################################################################################


$time=localtime();
use perlcyc;
$cyc = perlcyc -> new("ECOLI");

#Get arguments
$gene_links=$ARGV[0];
$dictionary=$ARGV[1];
$outfile=$ARGV[2];
chomp($outfile);

#Open OUT file and print header
open(OUT,">$outfile") || die "Cannot open OUT file at $outfile\n";
print OUT "# Created on $time from gene_reactions_flatfiles.pl\n# Fields:\n# [1] > Gene name\n# [2] > Reaction ID\n# [3] > Reaction direction\n# [4] > Substrates\n# [5] > Products\n";

###Get all genes
#Open gene links
open(GLK,$gene_links) || die "Cannot open GLK at $gene_links\n";
while(<GLK>){
  if($_=~/^([^\t]+)\t[^\t]+\t/){
    $gene_id=$1;
    
    #Get gene name
    $gene=&IDS($gene_id,"gi");

    #Get reactions
    @rxns=$cyc->reactions_of_gene($gene_id);
    foreach $rx (@rxns){

      #Reset variables
      $all_prods="";
      $all_subs="";

      #Get reaction direction
      @f_phyr=();
      @f_phyr=$cyc->get_slot_values($rx,"REACTION-DIRECTION");
      #If no reaction direction
      if(!(@f_phyr)){
	$f_phyr[0]="REVERSIBLE";
	print"Reaction $rx without direction. Considered reversible.\n";
      }

      #Find transport reactions
      $f_trans=0;
      $f_trans=$cyc->transport_rxn($rx);
      if($f_trans){
	if($f_phyr[0]=~/LEFT-TO-RIGHT/){
	  @subs=();
	  @subs=$cyc->get_slot_values($rx,"LEFT");
	  foreach $s (@subs){
	    if(!($s=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$|^\|*Pi\|*$/i)){
	      $s=&IDS($s,"sm");
	      $all_subs=$all_subs . "$s\_Ext//";
	    }
	  }
	  @prods=();
	  @prods=$cyc->get_slot_values($rx,"RIGHT");
	  foreach $p (@prods){
	    if(!($p=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$|^\|*Pi\|*$/i)){
	      $p=&IDS($p,"sm");
	      $all_prods=$all_prods . "$p//";
	    }
	  }
	  #Complete and save line
	  print OUT"$gene\t$rx\tL2R\t//$all_subs\t//$all_prods\n";
	}elsif($f_phyr[0]=~/RIGHT-TO-LEFT/){
	  @subs=();
	  @subs=$cyc->get_slot_values($rx,"RIGHT");
	  foreach $s (@subs){
	    if(!($s=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$|^\|*Pi\|*$/i)){
	      $s=&IDS($s,"sm");
	      $all_subs=$all_subs . "$s\_Ext//";
	    }
	  }
	  @prods=();
	  @prods=$cyc->get_slot_values($rx,"LEFT");
	  foreach $p (@prods){
	    if(!($p=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$|^\|*Pi\|*$/i)){
	      $p=&IDS($p,"sm");
	      $all_prods=$all_prods . "$p//";
	    }
	  }
	  print OUT"$gene\t$rx\tL2R\t//$all_subs\t//$all_prods\n";
	}elsif($f_phyr[0]=~/REVERSIBLE/){
	  @subs=();
	  @subs=$cyc->get_slot_values($rx,"LEFT");
	  foreach $s (@subs){
	    if(!($s=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$|^\|*Pi\|*$/i)){
	      $s=&IDS($s,"sm");
	      $all_subs=$all_subs . "$s\_Ext//";
	    }
	  }
	  @prods=();
	  @prods=$cyc->get_slot_values($rx,"RIGHT");
	  foreach $p (@prods){
	    if(!($p=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$|^\|*Pi\|*$/i)){
	      $p=&IDS($p,"sm");
	      $all_prods=$all_prods . "$p//";
	    }
	  }
	  print OUT"$gene\t$rx\tRVB\t//$all_subs\t//$all_prods\n";
	}
      }else{
	if($f_phyr[0]=~/LEFT-TO-RIGHT/){
	  @subs=();
	  @subs=$cyc->get_slot_values($rx,"LEFT");
	  foreach $s (@subs){
	    if(!($s=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$/i)){
	      $s=&IDS($s,"sm");
	      $all_subs=$all_subs . "$s//";
	    }
	  }
	  @prods=();
	  @prods=$cyc->get_slot_values($rx,"RIGHT");
	  foreach $p (@prods){
	  if(!($p=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$/i)){
	    $p=&IDS($p,"sm");
	    $all_prods=$all_prods . "$p//";
	  }
	}
	  print OUT"$gene\t$rx\tL2R\t//$all_subs\t//$all_prods\n";
	}elsif($f_phyr[0]=~/RIGHT-TO-LEFT/){
	  @subs=();
	  @subs=$cyc->get_slot_values($rx,"RIGHT");
	  foreach $s (@subs){
	    if(!($s=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$/i)){
	      $s=&IDS($s,"sm");
	      $all_subs=$all_subs . "$s//";
	    }
	  }
	  @prods=();
	  @prods=$cyc->get_slot_values($rx,"LEFT");
	  foreach $p (@prods){
	    if(!($p=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$/i)){
	      $p=&IDS($p,"sm");
	      $all_prods=$all_prods . "$p//";
	    }
	  }
	  print OUT"$gene\t$rx\tL2R\t//$all_subs\t//$all_prods\n";
	}elsif($f_phyr[0]=~/REVERSIBLE/){
	  @subs=();
	  @subs=$cyc->get_slot_values($rx,"LEFT");
	  foreach $s (@subs){
	    if(!($s=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$/i)){
	      $s=&IDS($s,"sm");
	      $all_subs=$all_subs . "$s//";
	    }
	  }
	  @prods=();
	  @prods=$cyc->get_slot_values($rx,"RIGHT");
	  foreach $p (@prods){
	    if(!($p=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$/i)){
	      $p=&IDS($p,"sm");
	      $all_prods=$all_prods . "$p//";
	    }
	  }
	  print OUT"$gene\t$rx\tRVB\t//$all_subs\t//$all_prods\n";
	}
      }
    }
  }
}
close(GLK);
close(OUT);


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

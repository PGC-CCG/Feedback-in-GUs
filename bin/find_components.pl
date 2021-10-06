################################################################################
##### 29-03-16                                                            ######
##### Find Components                                                     ######
##### Finds and quantifies components in get_short_paths_v3.pl output     ######
##### Improved from short_path_intersections.pl                           ######
##### Arguments:                                                          ######
##### [0]> short paths file without short path length restrictions        ######
##### [1]> Output file                                                    ######
##### Output:                                                             ######
##### [0] > TF                                                            ######
##### [1] > Total short paths                                             ######
##### [2] > Unique reactions from each component.                         ######
#####       Reactions separated by "," ; components separated by "//"     ######
##### History:                                                            ######
################################################################################

$time=localtime();

##Get args
$sp_file=$ARGV[0];
$outfile=$ARGV[1];
chomp($outfile);

undef %sps;

##CRetae temp file to keep alphabetical order in output
if($outfile=~/^(.+)\.txt/){
  $temp=$1 . "_temp.txt";
}else{
  die "Cannot create temp file, outfile is not .txt\n";
}

##Open out temp file
open(OUT,">$temp") || die "Cannot open OUT (temp) file at $temp.\n";
print OUT"#From find_components.pl on $time using:\n# $sp_file\n";

##Open IN file
open(IN,$sp_file) || die "Cannot open IN file at $sp_file.\n";
while(<IN>){
  if($_=~/^([^\t]+)\t([^\t]+)\t([^\t]+)\t/){
     $tf=$1;
     $id=$2 . "***" . $3;
     push(@{$sps{$tf}},$id);        
   }elsif($_=~/^([^\t]+)\t\t\t/){     # Fills index with all short paths per TF
     $tf=$1;
     print OUT"$tf\t0\n";
   }
}
close(IN);

#Create metabolite-short paths hash
foreach $tf (keys %sps){
  print"\n**$tf\n";
  undef %mets;
  $size1=(@{$sps{$tf}});      # Total short paths for that TF
  #Get metabolites
  for($i=0;$i<$size1;$i++){
    if($sps{$tf}[$i]=~/^(.+)\*\*\*(.+)/){
      $id=$1;           # Individual short path 
      @{$ids{$id}}=();
      $line=$2;
      @ssp=split(/\;\;/,$line);      # short path in metabolites: could be more than one
      foreach $ss (@ssp){    # for each possible short path
	@met=split(/\/\//,$ss);     # get individual metabolites
	foreach $m (@met){
	  if(!($m=~/^Pi$/)){
	    push(@{$mets{$m}},$id);      # {metabolite} (short path reactions) hash
	  }else{
	    print"$m was omitted from component search.\n";
	  }
	}
      }
    }else{
      print"$sps{$tf}[$i] was ommitted ; syntax in $sp_file is incompatible.\Find";
    }
  }
  
  #Find uniq reactions in short paths per metabolite 
  undef %uniq;
  foreach $met (keys %mets){
    @all_rxns=();
    foreach $ln (@{$mets{$met}}){
      @paths=split(/,/,$ln);
      @rxns=();
      foreach $shp (@paths){
	@rxns=split(/\/\//,$shp);
	push(@all_rxns,@rxns);
      }
    }

    $ref=&uniq(\@all_rxns);
    @all_rxns=@{$ref};
    push(@{$uniq{$met}},@all_rxns);
  }

  ### T E S T ###
  # Check array formation
  #foreach $ky (keys %uniq){
  #  print"$ky\t";
  #  foreach $el (@{$uniq{$ky}}){
  #    print"$el,";
  #  }
  #  print"\n";
  #}


  #Find groups
  @kys=keys %uniq;
  $size=@kys;
  foreach $ky1 (@kys){          # First key index
    $size1=@{$uniq{$ky1}};
 
   for($i=0;$i<$size1;$i++){    # First key array index
      $el1=$uniq{$ky1}[$i];
      
      for($j=0;$j<$size;$j++){    # 2nd key index
	$ky2=$kys[$j];
	if($ky1 eq $ky2){
	  next;
	}
	$size2=@{$uniq{$ky2}};

	for($q=0;$q<$size2;$q++){   #2nd key array index
	  $el2=$uniq{$ky2}[$q];

	  if($el1 eq $el2){   #if overlapping reactions among metabolites
	    push(@{$uniq{$ky1}},@{$uniq{$ky2}});      #join arrays to create group
	    @{$uniq{$ky2}}=();      #blank second array to avoid infinite loops
	    $ref=&uniq(\@{$uniq{$ky1}});               #delete repeated entries from joint array
	    @{$uniq{$ky1}}=@{$ref};                   #assign function value to array
	    $size1=@{$uniq{$ky1}};   #set first array new size
	    $q=$size2;               #break $q cycle and move to next key with the same 1st key array element


	    ### T E S T ###
	    #print"\n={$ky1}$el1 - {$ky2}$el2\n";
	    #foreach $ky (keys %uniq){
	    #  print"$ky\t";
	    #  foreach $el (@{$uniq{$ky}}){
	    #print"$el,";
	    #  }
	    #  print"\n";
	    #}

	  }
	}
      }
    }
  }

  #Print groups - print all keys with elements, each key should be a different group.
  $t_groups=0;
  foreach $ky (keys %uniq){
    $size=@{$uniq{$ky}};
    if($size){
      $t_groups++;
    }
  }
  print OUT"$tf\t$t_groups\t";
  foreach $ky (keys %uniq){
    $size=@{$uniq{$ky}};
    if($size){
      foreach $el (@{$uniq{$ky}}){
	print OUT"$el,";
      }
      print OUT"**";
    }
  }
  print OUT"\n";

}

close(OUT);

##Open out file and print header
open(LST,">$outfile") || die "Cannot open LST file at $outfile.\n";
print LST"#From find_components.pl on $time using:\n# $sp_file\n";
close(LST);

##Sort temp file
system("sort $temp >> $outfile");
system("rm $temp");




###############################
### subroutine uniq
### Returns uniq values from an array
### input: array
### output: uniq elements
###############################


sub uniq {

  local(@input)=@{$_[0]};
  local(@output)=();
  local($el,$new,$flag)=("","",0);

  foreach $el (@input){
    $flag=0;
    foreach $new (@output){
      if($el eq $new){
	$flag=1;
	last;
      }
    }
    if(!($flag)){
      push(@output,$el);
    }
  }
  
  return(\@output);

}



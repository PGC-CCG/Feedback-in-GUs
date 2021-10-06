#################################################################################
###### 20-11-20                                                            ######
###### Add object IDs                                                      ######
###### Adds object IDs to objects in a GU and propagates them to complexes.txt ##
###### modifications.txt, reactants_products.txt and objects.txt           ######
###### Argument usage:                                                     ######
###### [0] - directory with GUs path                                       ######
###### [1] - output directory path                                         ######
###### History:                                                            ######
###### Notes:                                                              ######
###### - Creates object ID with syntax obj+number                          ######
#################################################################################

use strict;
use warnings;
my $time=localtime();

###Get args
my $in_dir=$ARGV[0];
my $out_dir=$ARGV[1];
chomp($out_dir);

### Get all GU names
opendir(DIR,$in_dir) || die "Cannot open $in_dir.\n";
my @GUs=readdir(DIR);

foreach my $gu (sort @GUs){
  if($gu=~/^\./){
    next;
  }

  print"\n**$gu\n";
  #Get path name of GU files
  my $f_objects=$in_dir . "/" . $gu . "/" . "objects.txt";
  my $f_complexes=$in_dir . "/" . $gu . "/" . "complexes.txt";
  my $f_recprod=$in_dir . "/" . $gu . "/" . "reactants_products.txt";
  my $f_modification=$in_dir . "/" . $gu . "/" . "modification.txt";
  my $f_reactions=$in_dir . "/" . $gu . "/" . "reactions.txt";

  #Get path name of output GU files
  my $o_objects=$out_dir . "/" . $gu . "/" . "objects.txt";
  my $o_complexes=$out_dir . "/" . $gu . "/" . "complexes.txt";
  my $o_recprod=$out_dir . "/" . $gu . "/" . "reactants_products.txt";
  my $o_modification=$out_dir . "/" . $gu . "/" . "modification.txt";
  my $o_reactions=$out_dir . "/" . $gu . "/" . "reactions.txt";

  print"$o_complexes\n";
  
  #Create out folder and copy reactions.txt file
  my $new_dir=$out_dir . "/" . $gu;
  system("mkdir $new_dir");
  system("cp $f_reactions $o_reactions"); 
  
  ###Number objects and edit objects.txt
  # Create out objects file
  open(OOB,">$o_objects") || die "Cannot open $o_objects.\n";
  
  open(IN,$f_objects) || die "Cannot open IN at $f_objects.\n";
  my %objs_index;
  undef %objs_index;
  my $obid=0;
  while(<IN>){
    if($_=~/([^\t]+)\t([^\t]+)\t*([^\t]*)/){
      my $name=$1;
      my $type=$2;
      my $extra=$3;
      chomp($type);
      chomp($extra);
      $obid++;

      #Detect complexes and save ID
      if($type eq "COMPLEX"){
	#Find complex ID
	open(CPX,$f_complexes) || die "Cannot open CPX at $f_complexes.\n";
	while(<CPX>){
	  if($_=~/^(csa\d+)\t$name\t/){
	    $objs_index{$name}=$1;
	    print OOB"$1\t$name\t$type\t$extra\n";
	  }
	}
	close(CPX);
	$obid--;
	next;
      }
      my $key=$name;
      $objs_index{$key}="obj" . $obid;
      if($type=~/GENE|RNA/){
	print OOB"obj$obid\t$name\t$type\n";
      }else{
	print OOB"obj$obid\t$name\t$type\t$extra\n";
      }
    }
  }
  close(IN);
  close(OOB);

  ###Create new files
  #Rewrite complexes.txt
  open(CPX1,$f_complexes) || die "Cannot open $f_complexes.\n";
  open(OCX1,">$o_complexes") || die "Cannot open $o_complexes.\n";

  while(<CPX1>){
    if($_=~/^(csa\d+\t[^\t]+)\t([^\t]+)$/){
      my $item=$2;
      my $line=$1;
      chomp($item);

      print OCX1"$line\t$objs_index{$item}\t$item\n";
    }
  }
  close(CPX1);
  close(OCX1);

  #Rewrite modifications.txt
  open(MOD,$f_modification) || die "Cannot open $f_modification.\n";
  open(OMD,">$o_modification") || die "Cannot open $o_modification.\n";

  while(<MOD>){
    if($_=~/^(re\d+\t[^\t]+)\t([^\t]+)(\t*[^\t]*)$/){
      my $item=$2;
      my $line=$1;
      my $extra=$3;
      chomp($item);
      chomp($extra);

      if($extra=~/^\w/){
	print OMD"$line\t$objs_index{$item}\t$item\t$extra\n";
      }else{
	print OMD"$line\t$objs_index{$item}\t$item\n";
      }
    }
  }
  close(MOD);
  close(OMD);

  #Rewrite reactants_products.txt
  open(RPD,$f_recprod) || die "Cannot open $f_recprod.\n";
  open(ORP,">$o_recprod") || die "Cannot open $o_recprod.\n";

  while(<RPD>){
    if($_=~/^(re\d+)\t([^\t]+)\t([^\t]+)$/){
      my $item=$2;
      my $line=$1;
      my $extra=$3;
      chomp($extra);

      print ORP"$line\t$objs_index{$item}\t$item\t$extra\n";

    }
  }
  close(RPD);
  close(ORP);
  
}
closedir(DIR);





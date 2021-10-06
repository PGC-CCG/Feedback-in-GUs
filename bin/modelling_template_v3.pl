#################################################################################
###### 31-10-17                                                            ######
###### Modelling format                                                    ######
###### Creates files in a format than can be turned into equations         ######
###### [0] - gu_assembly directory                                         ######
###### [1] - Output dir                                                    ######
###### Output:                                                             ######
###### [1] > Objects                                                       ######
###### [2] > Effect                                                        ######
###### [3] > Object that is affected                                       ######
###### More info on: /Users/Daniela/Documents/PhD/Colaboraciones/          ######
######               Savageau/modelado/template_modelado_v2.txt            ######
###### History:                                                            ######
###### 29/11/17 - Included transport reactions. Omitted last comma in      ######
######            complex formation reactions.                             ######
###### Notes:                                                              ######
#################################################################################

use Sys::Hostname;
$host = hostname;
$time=localtime();
print"Begin: $time\n";

#Get args
my $gudir=$ARGV[0];
my $outdir=$ARGV[1];
chomp($outdir);

# Initializa out files
my $outname;
my $outfile;

#Open dir 
opendir(DIR,$gudir) || die "Cannot open DIR at $gudir.\n";
my @GU_files=readdir(DIR);
foreach my $file (@GU_files){
  if(!($file=~/^\./)){
    
    print"\n**$file\n";

    #Create outfile
    my $outfile=$outdir . $file . "_mod.txt";
    # Print header
    open(OUT,">$outfile") || die "Cannot open OUT at $outfile.\n";
    #print OUT"# From $0 on $time at $host\n# Input = $gudir\n";

    #Get file names from GU folder
    my $products=$gudir . $file . "/reactants_products.txt";
    my $objects=$gudir . $file . "/objects.txt";
    my $reactions=$gudir . $file . "/reactions.txt";
    my $modifications=$gudir . $file ."/modification.txt";
    my $complexes=$gudir . $file . "/complexes.txt";

    # Find catalysis reactions
    open(MOD,$modifications) || die "Cannot open MOD at $modifications.\n";
    while(<MOD>){
      if($_=~/^(re\d+)\tCATALYSIS\t([^\t]+)/){
	my $rex=$1;
	my $cat=$2;
	chomp($cat);

	my $f_transport=0;

	# Initialize variables for substrates & products
	my @substrates=();
	my @products=();

	# Get substrates and products
	open(RPD,$products) || die "Cannot open RPD at $products.\n";
	while(<RPD>){
	  if($_=~/^$rex\t([^\t]+)\treactant/){
	    my $rec=$1;
	    push(@substrates,$rec);
	  }elsif($_=~/^$rex\t([^\t]+)\tproduct/){
	    my $prod=$1;
	    push(@products,$prod);
	  }
	}
	close(RPD);
	
	# Find reversibility
	my $reverse="";
	open(REX,$reactions) || die "Cannot open REX at $reactions.\n";
	while(<REX>){
	  if($_=~/^$rex\tSTATE_TRANSITION\t([^\t]+)/){
	    $reverse=$1;
	    $f_transport=0;
	    last;
	  }
	  if($_=~/^$rex\tTRANSPORT\t([^\t]+)/){
	    $reverse=$1;
	    $f_transport=1;
	    last;
	  }
	}
	close(REX);


	# Print reaction

	if($reverse=~/^L2R$/){
	  # Print effect on all substrates
	  foreach my $rec1 (@substrates){
	    foreach my $rec2 (@substrates){
	      print OUT"$rec2,";
	    }
	    print OUT"$cat\t-\t$rec1";
	    if($f_transport){ 	# print type of reaction
	      print OUT"\tRT\n";
	    }else{
	      print OUT"\tRX\n";
	    }
	  }
	  # Print effect on all products
	  foreach my $prod1 (@products){
	    foreach my $rec2 (@substrates){
	      print OUT"$rec2,";
	    }
	    print OUT"$cat\t+\t$prod1";
	    if($f_transport){ 	# print type of reaction
	      print OUT"\tRT\n";
	    }else{
	      print OUT"\tRX\n";
	    }
	  }
	}elsif($reverse=~/^RVB$/){
	  # Print effect on all substrates on first direction
	  foreach my $rec1 (@substrates){
	    foreach my $rec2 (@substrates){
	      print OUT"$rec2,";
	    }
	    print OUT"$cat\t-\t$rec1";
	    if($f_transport){ 	# print type of reaction
	      print OUT"\tRTV\n";
	    }else{
	      print OUT"\tRXV\n";
	    }
	  }
	  # Print effect on all products on first direction
	  foreach my $prod1 (@products){
	    foreach my $rec2 (@substrates){
	      print OUT"$rec2,";
	    }
	    print OUT"$cat\t+\t$prod1";
	    if($f_transport){ 	# print type of reaction
	      print OUT"\tRTV\n";
	    }else{
	      print OUT"\tRXV\n";
	    }
	  }
	  # Print effect on all substrates on second direction
	  foreach my $rec1 (@products){
	    foreach my $rec2 (@products){
	      print OUT"$rec2,";
	    }
	    print OUT"$cat\t-\t$rec1";
	    if($f_transport){ 	# print type of reaction
	      print OUT"\tRTV\n";
	    }else{
	      print OUT"\tRXV\n";
	    }
	  }
	  # Print effect on all products on second direction
	  foreach my $prod1 (@substrates){
	    foreach my $rec2 (@products){
	      print OUT"$rec2,";
	    }
	    print OUT"$cat\t+\t$prod1";
	    if($f_transport){ 	# print type of reaction
	      print OUT"\tRTV\n";
	    }else{
	      print OUT"\tRXV\n";
	    }
	  }
	}
      }
    }
    close(MOD);

    # Print complexes
    my $last_cpx=1;
    my $i=0;
    # Find last complex
    open(CPX,$complexes) || die "Cannot open CPX at $complexes.\n";
    while(<CPX>){
      if($_=~/^csa(\d+)\t/){
	$last_cpx=$1;
      }
    }
    close(CPX);
    # Print all lines for each complex.
    for($i=1;$i<=$last_cpx;$i++){
      my @subunits=();
      my $csa_name="";
      my $monomer="";
      open(CPX,$complexes) || die "Cannot open CPX at $complexes.\n";
      while(<CPX>){
	if($_=~/^csa$i\t([^\t]+)\t([^\t]+)/){
	  $csa_name=$1;
	  $monomer=$2;
	  chomp($monomer);
	  push(@subunits,$monomer);
	}
      }
      close(CPX);
      # Print reactions type RC
      foreach my $sub1 (@subunits){
	my $psub=""; 		# Needed to eliminate commas at the EoL.
	foreach my $sub2 (@subunits){
	  $psub= $psub . "$sub2,";
	}
	chop($psub); 		# Eliminates last comma
	print OUT"$psub\t-\t$sub1\tRC\n";
      }
      my $psub2="";
      foreach my $sub2 (@subunits){
	$psub2= $psub2 . "$sub2,";
      }
      chop($psub2);
      print OUT"$psub2\t+\t$csa_name\tRC\n";
      # Print reactions type RCV
      foreach my $sub1 (@subunits){
	print OUT"$csa_name\t+\t$sub1\tRCV\n";
      }
    }

    # Find translation reactions
    open(RXN,$reactions) || die "Cannot open RXN at $reactions.\n";
    while(<RXN>){
      if($_=~/^(re\d+)\tTRANSLATION/){
	my $rex=$1;

	# Get TU and protein
	my $reactant="";
	my $product="";
	open(RPD,$products) || die "Cannot open RPD at $products\n";
	while(<RPD>){
	  if($_=~/^$rex\t([^\t]+)\treactant/){
	    $reactant=$1;
	  }elsif($_=~/^$rex\t([^\t]+)\tproduct/){
	    $product=$1;
	  }
	}
	close(RPD);

	# Print reaction
	print OUT"$reactant\t+\t$product\tTRL\n";
      }
    }
    close(RXN);

    # Find transcription reactions
    open(RXN,$reactions) || die "Cannot open RXN at $reactions.\n";
    while(<RXN>){
      if($_=~/^(re\d+)\tTRANSCRIPTION/){
	my $rex=$1;
	
	# Get transcribed mRNA
	my $product="";
	open(RPD,$products) || die "Cannot open RPD at $products\n";
	while(<RPD>){
	  if($_=~/^$rex\t([^\t]+)\tproduct/){
	    $product=$1;
	  }
	}
	close(RPD);

	# Get effect
	my $effect="";
	my $act_conf="";
	open(MOD,$modifications) || die "Cannot open MOD at $modifications.\n";
	while(<MOD>){
	  if($_=~/^$rex\t([^\t]+)\t([^\t]+)/){
	    $effect=$1;
	    $act_conf=$2;
	    chomp($act_conf);
	    # Print reaction
	    if($effect=~/^INHIBITION$/){
	      print OUT"$act_conf\t-\t$product\tTRC\n";
	    }elsif($effect=~/^PHYSICAL_STIMULATION$/){
	      print OUT"$act_conf\t+\t$product\tTRC\n";
	    }
	  }
	}
	close(MOD);

      }
    }
    close(RXN);
    close(OUT);
  }
}
closedir(DIR);


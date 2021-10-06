#################################################################################
###### 12-09-2020                                                          ######
###### Create compound dictionary file                                     ######
###### Retrieves protein, coding gene, common-name and synonyms            ######
###### Argument usage:                                                     ######
###### [0] - gene_links file                                               ######
###### [1] - gu_library path                                               ######
###### Output (depeding on type of molecule):                              ######
######  1)                                                                 ######
######   [0] > ProteinID                                                   ######
######   [1] > coding gene                                                 ######
######   [2] > //common-name//syn1//syn2//...//synX//                      ######
######  2)                                                                 ######
######   [0] > CompoundID                                                  ######
######   [1] > Common name                                                 ######
######   [2] > //syn1//syn2//...//synX//                                   ######
###### Notes:                                                              ######
###### - Requires Perlcyc & Pathway tools installed                        ######
###### - Pathway Tools should be running in -api mode                      ######
#################################################################################

$time=localtime();
use perlcyc;
$cyc = perlcyc -> new("ECOLI");

#Get arguments
$out_dir=$ARGV[0];
chomp($out_dir);

if($out_dir){
    $cmpdir=$out_dir . "gu_compound_dictionary.txt";
    open(OUT,">$cmpdir") || die "Cannot open OUT file at $cmpdir.\n";
}

print OUT"# From $0 on $time\n";

@proteins = $cyc -> get_class_all_instances('|Proteins|');
@subproteins = $cyc -> get_class_all_subs('|Proteins|');
@compounds = $cyc -> get_class_all_instances('|Compounds|');
@subcompounds = $cyc -> get_class_all_subs('|Compounds|');
@rnas = $cyc -> get_class_all_instances('|RNAs|');
@protein_complexes = $cyc -> get_class_all_instances('|Protein‑Complexes|');
@subprotein_complexes = $cyc -> get_class_all_subs('|Protein‑Complexes|');

foreach $prot(@proteins){
    @genes = $cyc -> get_slot_values($prot,'NAMES');
    if(scalar(@genes) > 1){
        $name = @genes[1];
    }
    else{
        $name = @genes[0];
    }
    print OUT"$prot\t$name\t//",join("//",@genes),"//\n";
}

foreach $subprot(@subproteins){
    @genes = $cyc -> get_slot_values($subprot,'NAMES');
    if(scalar(@genes) > 1){
        $name = @genes[1];
    }
    else{
        $name = @genes[0];
    }
    print OUT"$subprot\t$name\t//",join("//",@genes),"//\n";
}

foreach $protcom(@protein_complexes){
    $name = $cyc -> get_slot_value($ptocom,'NAME');
    @genes = $cyc -> get_slot_values($protcom,'NAMES');
    print OUT"$protcom\t$name\t//",join("//",@genes),"//\n";
}

foreach $subprotcom(@subprotein_complexes){
    $name = $cyc -> get_slot_value($subptocom,'NAME');
    @names = $cyc -> get_slot_values($subprotcom,'NAMES');
    $kegg_id = &get_KEGG_ID($subprotcom);
    if($kegg_id eq ""){
        print OUT"$subprotcom\t$name\t//",join("//",@names),"//\n";
    }
    else{
        print OUT"$subprotcom\t$name\t//$kegg_id//",join("//",@names),"//\n";
    }
    
}

foreach $cmpd(@compounds){
    $name = $cyc -> get_slot_value($cmpd,'COMMON-NAME');
    @names = $cyc -> get_slot_values($cmpd,'NAMES');
    $kegg_id = &get_KEGG_ID($cmpd);
    if($kegg_id eq ""){
        print OUT"$cmpd\t$name\t//",join("//",@names),"//\n";
    }
    else{
        print OUT"$cmpd\t$name\t//$kegg_id//",join("//",@names),"//\n";
    }
    
}

foreach $cmpdclass(@subcompounds){
    $name = $cyc -> get_slot_value($cmpdclass,'COMMON-NAME');
    @names = $cyc -> get_slot_values($cmpdclass,'NAMES');
    $kegg_id = &get_KEGG_ID($cmpdclass);
    if($kegg_id eq ""){
        print OUT"$cmpdclass\t$name\t//",join("//",@names),"//\n";
    }
    else{
        print OUT"$cmpdclass\t$name\t//$kegg_id//",join("//",@names),"//\n";
    }
}

foreach $rna(@rnas){
    @names = $cyc -> get_slot_values($rna,'NAMES');
    $kegg_id = &get_KEGG_ID($rna);
    if($kegg_id eq ""){
        print OUT"$rna\t@names[0]\t//",join("//",@names),"//\n";
    }
    else{
        print OUT"$rna\t@names[0]\t//$kegg_id//",join("//",@names),"//\n";
    }
}

sub get_KEGG_ID{

  local($eco_met)=($_[0]);
  my @dblinks=$cyc->get_slot_values($eco_met,'DBLINKS'); #Find slot with links to other databases.
  foreach my $kegg (@dblinks){
      if($kegg=~/error/){ #Omit if slot is not found (returns an error message)
        print"$eco_met\tERROR\n";
        last;
      }

      if(@{$kegg}[0]=~/LIGAND-CPD/){ #Slot is an array with a database per item. This identifies the element in @kegg that contains the KEGG IDs
        my $kegg_id=@{$kegg}[1];
        return($kegg_id);
        #print OUT2"$eco_met\t@{$kegg}[0]\t@{$kegg}[1]\n"; #Prints the database identifier (should read LIGAND-CPD)
      }
    }

}

close(OUT);
#################################################################################
###### 07-20-2020                                                          ######
###### Author: Georgette Femerling                                         ######
######                                                                     ######
###### Gets compounds of pathways and gets pathways with a distance of 0   ######
###### Argument usage:                                                     ######
###### [0] - path to gu_library directory                                  ######
###### Output:                                                             ######
######  1 - pathway_compounds.txt :                                        ######
######      [1] Pathway                                                    ######
######      [2] compounds                                                  ######
######  2 - pathways_cerodistance.txt                                      ######
######      [1] pathway1                                                   ######
######      [2] pathway2                                                   ######
######      [3] shared compounds                                           ######
######      [4] # of shared compounds                                      ######
######  3 - pathway_enzymes.txt :                                          ######
######      [1] Pathway                                                    ######
######      [2] enzymes                                                    ######
###### History:                                                            ######
#################################################################################

#Load modules
use perlcyc;
$cyc = perlcyc -> new("ECOLI");

#get arguments
$time=localtime();
$dir_out=$ARGV[0];
chomp($dir_out);

#create output files
$output=$dir_out . "pathway_compounds.txt";
$output2=$dir_out . "pathways_cerodistance.txt";
$output3=$dir_out . "pathway_enzymes.txt";

#Open outfiles
open(OUT,">$output") || die "Cannot open OUT file at $output.\n";
open(OUT1,">$output2") || die "Cannot open OUT1 file at $output2.\n";
open(OUT2,">$output3") || die "Cannot open OUT2 file at $output3.\n";

#Print output files headers
print OUT"# From $0 on $time\n# Pathway\tCompounds\n";
print OUT1"# From $0 on $time\n# Pathway1\tPathway2\tShared compounds\tNumbber of shared compounds\n";
print OUT2"# From $0 on $time\n# Pathway\tEnzymes\n";

#Get all ecoli pathways 
@pathways = $cyc -> all_pathways();

my %pwys;

#Get compounds of each pathway

foreach $pwy (@pathways){
    #print OUT"$pwy\t";
    @compounds = $cyc -> compounds_of_pathway($pwy);
    @{$pwys{$pwy}} = @compounds;
    print OUT "$pwy\t";
    print OUT join(",", @{$pwys{$pwy}}),"\n";
    @enzymes = $cyc -> enzymes_of_pathway($pwy);
    print OUT2 "$pwy\t";
    print OUT2 join(",", @enzymes),"\n";
}

#Generate combinations between pathways
@pathways2 = @pathways;

foreach $pwy (@pathways){
    foreach $pwy2 (@pathways2){
        my $count = 0;
        my @intersect;
        if($pwy eq $pwy2){
            next;
        }
        foreach $element (@{$pwys{$pwy}}){
            if($element=~/^PROTON$|^.MP$|^.TP$|^.DP$|^WATER$|^NAD.*$|^\|*Pi\|*$|^OXYGEN-MOLECULE$|^CARBON-DIOXIDE$|^CO-A$|^PPI$/i){
                next;
            }
            if (grep { $_ eq $element} @{$pwys{$pwy2}}){
                push(@intersect, $element);
                $count++;
            }
        }
        if(@intersect){
            print OUT1"$pwy\t$pwy2\t";
            print OUT1 join(",",@intersect);
            print OUT1"\t$count\n";
        }
        #Eliminar el pwy que ya se comparo con todo del segundo 
        @pathways2 = grep {$_ ne $pwy} @pathways2;
    }
}

close(OUT);
close(OUT1);
close(OUT2);

print "----------- DONE ------------\n";
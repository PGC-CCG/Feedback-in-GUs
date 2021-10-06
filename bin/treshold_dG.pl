use perlcyc;
$cyc = perlcyc -> new("ECOLI");

@all_rxs = $cyc -> all_rxns();

foreach $rx(@all_rxs){
    $dg = $cyc -> get_slot_value($rx, "GIBBS-0");
    $dir = $cyc -> get_slot_value($rx,"REACTION-DIRECTION");
    print $dg,"\t",$dir,"\n";
}
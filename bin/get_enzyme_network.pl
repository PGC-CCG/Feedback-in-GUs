
#Load Modules
#export PERL5LIB=/home/gfemer/GUs/
use perlcyc;
$cyc = perlcyc -> new("ECOLI");

#Get all enzymes

@enzymes = $cyc -> all_enzymes();
#print join("\n",@enzymes);

#Get reactions of enzymes
foreach $enzyme (@enzymes){
    #print "------------$enzyme-----------\n";
    @rxs = $cyc -> reactions_of_enzyme($enzyme);
    foreach $rx (@rxs){
        #print "$enzyme\t$rx\t";
        #@substrates = $cyc -> substrates_of_reaction($rx);
        #get direction of reaction
        $direction = $cyc -> get_slot_value($rx,"REACTION-DIRECTION");
        #print "$direction\t";
        if ($direction=~/RIGHT-TO-LEFT/i){
            print "$enzyme\t";
            @subs=$cyc->get_slot_values($rx,"LEFT");
            @prods=$cyc->get_slot_values($rx,"RIGHT");
            print join(",",@subs),"\t";
            print join(",",@prods),"\n";
        }
        elsif ($direction=~/LEFT-TO-RIGHT/i){
            $direc = "L2R";
            print "$enzyme\t";
            @subs=$cyc->get_slot_values($rx,"RIGHT");
            @prods=$cyc->get_slot_values($rx,"LEFT");
            print join(",",@subs),"\t";
            print join(",",@prods),"\n";
        }
        elsif ($direction=~/REVERSIBLE/i){
            print "$enzyme\t";
            @subs=$cyc->get_slot_values($rx,"LEFT");
            @prods=$cyc->get_slot_values($rx,"RIGHT");
            print join(",",@subs),",",join(",",@prods),"\t";
            print join(",",@subs),",",join(",",@prods),"\n";
        }
    }
}
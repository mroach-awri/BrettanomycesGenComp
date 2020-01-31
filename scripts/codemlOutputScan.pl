use strict; 
use warnings;

# one-off script to return BEB significant hits

my $fofn = shift or die;

open my $FO, '<', $fofn or die;

while(my $f = <$FO>){
    chomp($f);
    my @o;
    my $score = `cat $f`;
    $score =~ s/^[\s\S]*BEB[\s\S]*SE for w//;
    #$score =~ s/^[\s\S]*NEB[\s\S]*SE for w//;
    $score =~ s/The grid[\s\S]*$//;
    #$score =~ s/Bayes Empirical[\s\S]*$//;
    my@l=split/\n/,$score;
    for my $l (@l){
        if ($l =~ m/\*/){
            $l =~ s/^\s+//;
            $l =~ s/\s/\t/;
            push @o, "$l\n";
        }
    }
    
    if (@o > 0){
        $f =~ s/.*OG/OG/;
        $f =~ s/\..*//;
        print "$f\n@o";
    }
}

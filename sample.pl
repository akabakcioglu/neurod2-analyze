#/usr/bin/perl -w
use strict;

my $n = shift @ARGV;		# number of samples to be picked
my @samples = ();

my $line;
foreach my $i (1..$n) {
    $line = <> or die "Not enough lines in the input.\n";
    push @samples,$line 
}

my $cnt = $n;
while ($line = <>) {
    $cnt += 1;
    if (rand() < (1.0*$n)/$cnt) {
	my $i = int(rand($n));
	$samples[$i] = $line;
    }
}

foreach my $l (@samples) {
    print $l;
}

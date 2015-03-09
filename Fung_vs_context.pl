#/usr/bin/perl -w
use strict;
my $Funglist = $ARGV[0];
open(FP,$Funglist) or die "Can't open $Funglist\n";
my %h=();
while(<FP>) {
    chomp;
    $h{$_}=1;
}
close(FP);
print STDERR "number of proteins = " . scalar(keys %h) . "\n";
my $ourlist = $ARGV[1];
open(FP,$ourlist) or die "Can't open $ourlist\n";
my $n=0;
while(<FP>) {
    chomp;
    my $p = $_;
#    my @a = split /,/;
#    foreach my $x (@a) {
#    if (defined $h{lc($x)}) {
    if (defined $h{$p}) {
	$n += 1;
#	}
    }
}
close(FP);
print STDERR "overlap = $n.\n";



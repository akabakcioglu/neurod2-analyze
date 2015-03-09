#/usr/bin/perl -w
use strict;
use List::MoreUtils qw(uniq);

$|=1;				# flush print commands

my $fname = shift @ARGV;
open(FP,$fname) or die "Run as: > perl $0.pl PATH/bindingregion_*.txt > outfile.\n";

my %h=();
while (<FP>){			# read one sequence at a time
    chomp;
    my ($reg,$ord,$nhits,$len_avg,$zscore,$peakscore_avg,$genelist,$peaklist) = split/\t/;
    foreach my $x (split/,/,$peaklist) {
	$h{$x}{$reg}+=1;
    }
}
close(FP);

my %hh=();			# get venn diagram entries
foreach my $x (keys %h) {
    push @{$hh{join(".",sort keys %{$h{$x}})}},$x;
}

foreach my $x (sort keys %hh) {
    print STDOUT join("\t",$x,scalar @{$hh{$x}},join(",",sort @{$hh{$x}}))."\n";
}


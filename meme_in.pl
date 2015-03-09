#/usr/bin/perl -w
use strict;

$|=1;				# flush print commands

my $reg = shift @ARGV;
my $order = shift @ARGV;
my $f1 = shift @ARGV;
my $f2 = shift @ARGV;

die "Run as: > perl $0.pl (utr3|utr5|cds|intron) (0|1|..) (bindingregion_*.txt) (merge_seq.txt) > outfile.\n" 
    unless ($reg =~ /utr3|utr5|cds|intron/ && 
	    $order =~ /^[0-9]+$/);

print STDERR "Reading peak list... ";
open(FP,$f1) or die "Can not open $f1.\n";
my $line;
do {
    $line = <FP>;
} while ($line !~ /^$reg\t$order\t/);
chomp($line);
my @a = split /\t/,$line;
my @seqlist = split /,/,$a[-1];
close(FP);
my $listsize = scalar(@seqlist);
if ($listsize==0) {
    die "No peaks found peaks.\n";
} else {
    print STDERR "got $listsize peaks.\n";
}

my %h=();
foreach my $x (@seqlist) {
    $h{$x}=1;
}

print STDERR "Reading sequences... ";
open(FP,$f2) or die "Can not open $f2.\n";
while(<FP>){
    if (/^>([0-9]*)$/) {
	my $seq = <FP>;
	print ">$1\n$seq" if (defined $h{$1});
    } else {
	die "merge_seq.txt file not as expected.\n";
    }
}



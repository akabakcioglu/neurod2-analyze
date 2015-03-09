#/usr/bin/perl -w
use strict;

$|=1;				# flush print commands

my %id2gene=();
$id2gene{"N/A"}="N/A";
my $filename = "/home/alkan/c/gulayse/neurod2/mouse_data/ucsc/ensemblToGeneName.txt";
open(FP,$filename) or die "Can't open ensemblToGeneName.txt";
while(<FP>) {
    chomp;
    my ($id,$genename) = split/\t/;
    $id2gene{$id}=$genename;
}


my $mousepath = "/home/alkan/c/gulayse/neurod2/mouse_data/ucsc/";
my $filename = "ensGene.txt";
open(FP,$mousepath.$filename) or die "Can't open $filename";

my %hgenelen=();

while(<FP>) {
    chomp;
    my ($chr,$txstart,$txend,$id,$a,$sense,@unused) = split/\t/;
    my $l = $txend-$txstart;
    unless (defined $hgenelen{$id2gene{$id}} && $hgenelen{$id2gene{$id}} > $l) {
	$hgenelen{$id2gene{$id}} = $l; # pick max length tx for this gene name
    }
}
close(FP);
print STDERR "finished reading $filename.\n";


my %hscore=();
my %hitcount=();

my $margin = 10000;
my $fmargin = 10000.0;

my $tightmargin = 1000;


my $cnt = 0;
while (<>){			# read context_ens.txt
    $cnt += 1;
    print STDERR "." if ($cnt % 100 == 0);
    chomp;
    my ($id,$score,$ingene,$inlen,$uptssgene,$uptsslen,$downtssgene,$downtsslen,$downgene,$downlen,$ensliststring) = split/\t/;
    my @enslist = split/,/,$ensliststring;

    if ($ingene !~ /N\/A/ && $inlen =~ /[0-9]*/) {
	$hscore{$ingene}{in} += $fmargin/($fmargin+$hgenelen{$ingene});
	# $hscore{$ingene}{in} += 1
	#     if ($downtsslen <= $tightmargin);
	$hitcount{$ingene} += 1;
    }
    if (-$uptsslen <= $margin) {	# $uplen is negative
	$hscore{$uptssgene}{up} += 1; # unless (defined $h{$uptssgene} && $h{$uptssgene} > 1);
	# $hscore{$uptssgene}{up} += 1
	#     if (-$uptsslen <= $tightmargin);
	#     if (defined $hscore{$uptssgene}{down} and $hscore{$uptssgene}{down} > $hscore{$uptssgene}{up});
	$hitcount{$uptssgene} += 1;
    }
    if ($downlen <= $margin) {	# $uplen is negative
	$hscore{$downgene}{down} += 1; # unless (defined $h{$downgene} && $h{$downgene} > 1);
	# $hscore{$downgene}{down} += 1
	#     if ($downlen <= $tightmargin);
	#     if (defined $hscore{$downgene}{up} and $hscore{$downgene}{up} > $hscore{$downgene}{down});
	$hitcount{$downgene} += 1;
    }
}    
my %htotal = ();
for my $x (keys %hscore) {
    $htotal{$x} = $hscore{$x}{in} + $hscore{$x}{up} + $hscore{$x}{down};
#    $htotal{$x} *= 2 if ($hscore{$x}{up} > 0 && $hscore{$x}{down} > 0);
}

for my $x (reverse sort {$htotal{$a}<=>$htotal{$b}} keys %htotal) {
    print "$x\t$htotal{$x}\t$hgenelen{$x}\t$hitcount{$x}\t$hscore{$x}{in}\t$hscore{$x}{up}\t$hscore{$x}{down}\n" if ($hitcount{$x} > 1 && $hgenelen{$x}>1000);
}

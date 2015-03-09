#/usr/bin/perl -w
use strict;

$|=1;				# flush print commands

my $VERYLARGENUMBER=100000000;
my %hpeak=();
my %hchr=();

my $filetype = shift @ARGV;
my $offset = shift @ARGV;

die "Run as: > perl $0.pl (ens|refseq) (offset) < peakfile > outfile.\n" 
    unless ($filetype =~ /ens|refseq/);
print STDERR "Reading peak data file... ";
while (<>){			# read one sequence at a time
    chomp;
    my ($id,$nmerge,$score,$chr,$start,$end,$lmerge,$mstart,$mend,@mergelist) = split/\t/;
    $hpeak{$id}{score} = $score;
    $hpeak{$id}{chr} = $chr;
    $hpeak{$id}{center} = int(0.5*($start+$end));
    $hchr{$chr}{$id}=1;
}
my $npeak = keys %hpeak;
print STDERR "finished reading $npeak peaks.\n";

# Read each chromosome in mm10 and locate the peaks within genes
my $pathname = "/home/alkan/c/gulayse/neurod2/mouse_data/ucsc/";
my $filename;

if ($filetype =~ /ens/) {
    $filename = "ensGene.txt";
} elsif ($filetype =~ /refseq/) {
    $filename = "RefSeq_refgene.txt";
}
open(FP,$pathname.$filename) or die "Can't open $filename.\n";

print STDERR "Reading transcript data ($filename) and locating peaks ..\n";
my %hgene = ();
## remove the header if there is one
<FP>; unless (/^\#/) {close(FP); open(FP,$pathname.$filename) or die "Can't open $filename.\n";} # reset pointer

my @unused=();
my ($chr,$txstart,$txend,$id,$sense,$bin,$cdsstart,$cdsend,$exoncnt,$exonstarts,$exonends,$b,$genename);

my $startcnt=0;
my $endcnt=0;
my $ctrcnt=0;
my $startend = 0;
my $startctr = 0;
my $endctr = 0;

while(<FP>) {
    chomp;
    if ($filetype =~ /ens/) {	# checking again and again is stupid, but let's not worry about it
	($chr,$txstart,$txend,$id,$a,$sense,@unused) = split/\t/;
    } elsif ($filetype =~ /refseq/) {
	($bin,$id,$chr,$sense,$txstart,$txend,$cdsstart,$cdsend,$exoncnt,$exonstarts,$exonends,$b,$genename,@unused) = split/\t/;
    }

    my $hitstart=0;
    my $hitend = 0;
    my $hitctr = 0;
    my $txctr = 0.5*($txend+$txstart);

    for my $pid (keys %{$hchr{$chr}}) { # check if any peak falls here
	my $ctr = $hpeak{$pid}{center};
	if ($ctr > $txstart-$offset && $ctr < $txstart+$offset) { # peak center near tx start
	    if ($sense eq "+") {
		$hitstart += 1;
		$startcnt += 1;
	    } elsif ($sense eq "-") {
		$hitend += 1;
		$endcnt += 1;
	    }
	} 
	if ($ctr > $txend-$offset && $ctr < $txend+$offset) { # peak center near tx end
	    if ($sense eq "-") {
		$hitstart += 1;
		$startcnt += 1;
	    } elsif ($sense eq "+") {
		$hitend += 1;
		$endcnt += 1;
	    }
	}
	if ($ctr > $txctr-$offset && $ctr < $txctr+$offset) { # peak center near tx center
	    $hitctr += 1;
	    $ctrcnt += 1;
	}
    }

    if ($hitstart && $hitend) {	# hit both ends
	$startend += 1;
    } elsif ($hitstart && $hitctr) {
	$startctr += 1;
    } elsif ($hitend && $hitctr) {
	$endctr += 1;
    }
}

print STDERR "start-hits: $startcnt\n";
print STDERR "end-hits: $endcnt\n";
print STDERR "center-hits: $ctrcnt\n";
print STDERR "start+end-hits: $startend\n";
print STDERR "start+center-hits: $startctr\n";
print STDERR "end+center-hits: $endctr\n";
    
close(FP);


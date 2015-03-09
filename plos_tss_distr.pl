#/usr/bin/perl -w
use strict;

$|=1;				# flush print commands

my $bound=1000000;
my %hpeak=();
my %hchr=();

my @updist=();
my @downdist=();
foreach my $i (0..$bound) {$updist[$i]=0; $downdist[$i]=0;}
my $upnorm = 0;
my $downnorm = 0;

my $filetype = shift @ARGV;
die "Run as: > perl $0.pl (ens|refseq) < inputfile(peaks) > outfile(cummulative_prob_distr). Specified type must be either ens or refseq!\n" 
    unless ($filetype =~ /ens|refseq/);
print STDERR "Reading peak data file... ";
while (<>){			# read one sequence at a time
    chomp;
    my ($id,$nmerge,$score,$chr,$start,$end,$lmerge,$mstart,$mend,@mergelist) = split/\t/;
    $hpeak{$id}{score} = $score;
    $hpeak{$id}{chr} = $chr;
    $hpeak{$id}{start} = $start;
    $hpeak{$id}{end} = $end;
    $hpeak{$id}{center} = int(0.5*($start+$end));
    $hchr{$chr}{$id}=1;
}
my $npeak = keys %hpeak;
print STDERR "finished reading $npeak peaks.\n";

# Read each chromosome in mm10 and locate the peaks within genes
my $pathname = "/home/alkan/c/gulayse/neurod2/mouse_data/ucsc/";
my $filename;

if ($filetype =~ /ens/) {
    $filename = "ensGene_uniq.txt";
} elsif ($filetype =~ /refseq/) {
    $filename = "RefSeq_refgene_uniq.txt";
}

open(FP,$pathname.$filename) or die "Can't open $filename.\n";

print STDERR "Reading transcript data ($filename) and locating peaks ..\n";
## remove the header if there is one
<FP>; unless (/^\#/) {close(FP); open(FP,$pathname.$filename) or die "Can't open $filename.\n";} # reset pointer

my @unused=();
my ($chr,$txstart,$txend,$id,$sense,$bin,$cdsstart,$cdsend,$exoncnt,$exonstarts,$exonends,$a,$genename);

while(<FP>) {
    chomp;
    if ($filetype =~ /ens/) {	# checking again and again is stupid, but let's not worry about it
	($chr,$txstart,$txend,$id,$a,$sense,@unused) = split/\t/;
    } elsif ($filetype =~ /refseq/) {
	($bin,$id,$chr,$sense,$txstart,$txend,@unused) = split/\t/;
    }

    for my $pid (keys %{$hchr{$chr}}) { # check if any peak falls here
	my $ctr = $hpeak{$pid}{center};
	my $tssloc = $bound+1;
	if ($sense eq "+") {
	    $tssloc = $ctr-$txstart;
	} elsif ($sense eq "-") {
	    $tssloc = -$ctr+$txend;
	}
	my $upstream = $tssloc < 0 ? 1 : 0;
	my $tssdist = $tssloc >= 0 ? $tssloc : -$tssloc;
	unless ($tssloc < -$bound or $tssloc > $bound){ # within +/- 10^6 of tss
	    if ($upstream == 1) {
		$updist[$tssdist] += 1;
		$upnorm += 1;
	    } else {
		$downdist[$tssdist] += 1;
		$downnorm += 1;
	    }
	}
    }				# end of for-loop over peaks
}
close(FP);

foreach my $i (0..$bound) {
    $updist[$i] /= $upnorm;
    $downdist[$i] /= $downnorm;
}

my $bin = 20;
my $cnt = 0;
while (@updist) {
    my $valup=0;
    my $valdown=0;
    foreach my $x (1..$bin) {
	$valup += shift @updist;
	$valdown += shift @downdist;
    }
    $cnt++;
    print join("\t",$cnt,$valup,$valdown)."\n";
}

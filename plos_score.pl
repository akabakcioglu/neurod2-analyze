#/usr/bin/perl -w
use strict;

$|=1;				# flush print commands

my $bound=1000000;
my %hpeak=();
my %hchr=();
my %id2gene=();

my %closest=();
my %closestdist=();
my %isup=();

my @updist=();
my @downdist=();
foreach my $i (0..$bound) {$updist[$i]=0; $downdist[$i]=0;}
my $upnorm = 0;
my $downnorm = 0;

my $filetype = shift @ARGV;
die "Run as: > perl $0.pl (ens|refseq) < inputfile(peaks) > outfile (tx scores). Specified type must be either ens or refseq!\n" 
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

if ($filetype =~ /ens/) {	# ens id to genename mapping is in a separate file
    my $filename = "/home/alkan/c/gulayse/neurod2/mouse_data/ucsc/ensemblToGeneName.txt";
    open(FP,$filename) or die "Can't open ensemblToGeneName.txt";
    while(<FP>) {
	chomp;
	my ($id,$geneid) = split/\t/;
	$id2gene{$id}=$geneid;
    }
    close(FP);
}

open(FP,$pathname.$filename) or die "Can't open $filename.\n";

print STDERR "Reading transcript data ($filename) and locating peaks ..\n";
## remove the header if there is one
<FP>; unless (/^\#/) {close(FP); open(FP,$pathname.$filename) or die "Can't open $filename.\n";} # reset pointer

my @unused=();
my ($chr,$txstart,$txend,$id,$sense,$bin,$cdsstart,$cdsend,$exoncnt,$exonstarts,$exonends,$ab,$genename);

while(<FP>) {
    chomp;
    if ($filetype =~ /ens/) {	# checking again and again is stupid, but let's not worry about it
	($chr,$txstart,$txend,$id,$ab,$sense,@unused) = split/\t/;
    } elsif ($filetype =~ /refseq/) {
	($bin,$id,$chr,$sense,$txstart,$txend,$cdsstart,$cdsend,$exoncnt,$exonstarts,$exonends,$ab,$genename,@unused) = split/\t/;	
	$id2gene{$id}=$genename;
    }
    next if ($id2gene{$id} =~ /^Gm|^Mir|Rik$|\-mir\-/);
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
	unless ($tssdist > $bound){ # within +/- 10^6 of tss
	    if ($upstream == 1) {
		$updist[$tssdist] += 1;
		$upnorm += 1;
	    } else {
		$downdist[$tssdist] += 1;
		$downnorm += 1;
	    }
	    unless (defined $closestdist{$pid} and $closestdist{$pid} < $tssdist) {
		$closestdist{$pid} = $tssdist;
		$closest{$pid} = $id;
		$isup{$pid} = $upstream;
	    }
	}
    }				# end of for-loop over peaks
}
close(FP);

$updist[0] /= $upnorm;
$downdist[0] /= $downnorm;
foreach my $i (1..$bound) {
    $updist[$i] = $updist[$i-1] + ($updist[$i]/$upnorm);
    $downdist[$i] = $downdist[$i-1] + ($downdist[$i]/$downnorm);
}

# now get transcript scores.

my %hscore=();
my %hpeaksof=();
my %hclosestof=();

my $npeak_isolated = 0;
for my $pid (keys %hpeak) {
    if (defined $isup{$pid}) { 
	if ($isup{$pid}==1) {
	    $hscore{$closest{$pid}} += -log($updist[$closestdist{$pid}]);
	    $hpeaksof{$closest{$pid}}{$pid}=1;
	    # record the nearest (upstream) peak to the gene and its distance
	    unless (defined $hclosestof{$closest{$pid}}{up}{dist}
		    and $hclosestof{$closest{$pid}}{up}{dist} < $closestdist{$pid}) {
		$hclosestof{$closest{$pid}}{up}{pid} = $pid;
		$hclosestof{$closest{$pid}}{up}{dist} = $closestdist{$pid};
	    }
	} elsif ($isup{$pid}==0) {
	    $hscore{$closest{$pid}} += -log($downdist[$closestdist{$pid}]);
	    $hpeaksof{$closest{$pid}}{$pid}=-1;
	    # record the nearest (downstream) peak to the gene and its distance
	    unless (defined $hclosestof{$closest{$pid}}{down}{dist}
		    and $hclosestof{$closest{$pid}}{down}{dist} < $closestdist{$pid}) {
		$hclosestof{$closest{$pid}}{down}{pid} = $pid;
		$hclosestof{$closest{$pid}}{down}{dist} = $closestdist{$pid};
	    }
	} else {
	    die "value other than 0,1 for $pid in \%isup !!\n";
	}
    } else {
	$npeak_isolated += 1;
    }
}
print STDERR "$npeak_isolated isolated peaks.\n";

for my $txid (reverse sort {$hscore{$a}<=>$hscore{$b}} keys %hscore) { # check if any peak falls here
    print join("\t",$txid,$hscore{$txid},$id2gene{$txid},$hclosestof{$txid}{up}{pid},$hclosestof{$txid}{up}{dist},$hclosestof{$txid}{down}{pid},$hclosestof{$txid}{down}{dist})."\n";
}

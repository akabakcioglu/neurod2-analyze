#/usr/bin/perl -w
use strict;

$|=1;				# flush print commands

my $VERYLARGENUMBER=100000000;
my %hpeak=();
my %hchr=();

my $filetype = shift @ARGV;
die "Run as: > perl $0.pl (ens|refseq) < inputfile > outfile. Specified type must be either ens or refseq!\n" 
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
    $hpeak{$id}{tssloc} = $VERYLARGENUMBER;
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
## remove the header if there is one
<FP>; unless (/^\#/) {close(FP); open(FP,$pathname.$filename) or die "Can't open $filename.\n";} # reset pointer

my @unused=();
my ($chr,$txstart,$txend,$id,$sense,$bin,$cdsstart,$cdsend,$exoncnt,$exonstarts,$exonends,$b,$genename);

while(<FP>) {
    chomp;
    if ($filetype =~ /ens/) {	# checking again and again is stupid, but let's not worry about it
	($chr,$txstart,$txend,$id,$a,$sense,@unused) = split/\t/;
    } elsif ($filetype =~ /refseq/) {
	($bin,$id,$chr,$sense,$txstart,$txend,$cdsstart,$cdsend,$exoncnt,$exonstarts,$exonends,$b,$genename,@unused) = split/\t/;
    }

    for my $pid (keys %{$hchr{$chr}}) { # check if any peak falls here
	my $ctr = $hpeak{$pid}{center};
	my $tssloc = $ctr-$txstart;
	if (abs($tssloc) < abs($hpeak{$pid}{tssloc})){
	    $hpeak{$pid}{tssloc} = $tssloc;
	    $hpeak{$pid}{tssloc} *= -1 if ($sense eq "-");
	}
    }				# end of for-loop over peaks
}
close(FP);
print STDERR "done extracting peak-tss distances. Now converting ensembl IDs to gene names..\n";

print STDERR "printing data: merge_id / score / location relative to tss..\n";

for my $pid (keys %hpeak) {
    print STDOUT join("\t",
		      $pid,
		      $hpeak{$pid}{score},
		      $hpeak{$pid}{tssloc}
	)."\n";
}
print STDERR "done.\n";


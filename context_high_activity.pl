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
    @{$hpeak{$id}{intx}} = ();
    $hpeak{$id}{upgene} = "N/A";
    $hpeak{$id}{downgene} = "N/A";
    $hpeak{$id}{uploc} = $VERYLARGENUMBER;
    $hpeak{$id}{downloc} = $VERYLARGENUMBER;
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

my %id2gene=();
$id2gene{"N/A"}="N/A";
if ($filetype =~ /ens/) {	# ens id to genename mapping is in a separate file
    my $filename = "/home/alkan/c/gulayse/neurod2/mouse_data/ucsc/ensemblToGeneName.txt";
    open(FP,$filename) or die "Can't open ensemblToGeneName.txt";
    while(<FP>) {
	chomp;
	my ($id,$geneid) = split/\t/;
	$id2gene{$id}=$geneid;
    }
}

my %genelen = ();
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
	$id2gene{$id}=$genename;
    }
    my $glen = $txend-$txstart;
    my $genename = $id2gene{$id};
    $genelen{$genename} = $glen unless (defined $genelen{$genename} && $genelen{$genename} >= $glen);

    for my $pid (keys %{$hchr{$chr}}) { # check if any peak falls here
##	my $st = $hpeak{$pid}{start};
##	my $en = $hpeak{$pid}{end};
##	if ($st > $txstart && $en < $txend) { # peak in transciption region
	my $ctr = $hpeak{$pid}{center};
	if ($ctr >= $txstart && $ctr <= $txend) { # peak center in transciption region
	    push @{$hpeak{$pid}{intx}},$id;
	}
	# check up/down stream of tx-start
	if ($sense eq "+") { # upstream is smaller than $txstart
	    if ($ctr < $txstart) { # peak is upstream of tss
		if ($hpeak{$pid}{uploc} > ($txstart-$hpeak{$pid}{center})) { # update "upstream of tss" info
		    $hpeak{$pid}{upgene} = $id;
		    $hpeak{$pid}{uploc} = $txstart-$hpeak{$pid}{center};
		    }
	    } elsif (($hpeak{$pid}{center}-$txend) > 0) {		# peak is downstream of gene-end
		if ($hpeak{$pid}{downloc} > ($hpeak{$pid}{center}-$txend)) { # update "downstream of gene-end" info
		    $hpeak{$pid}{downgene} = $id;
		    $hpeak{$pid}{downloc} = $hpeak{$pid}{center}-$txend;
		    }
	    }
	} else {	      # if sense is "-" switch "up" and "down"
	    if ($hpeak{$pid}{center} < $txstart) { # peak is downstream of tss
		if ($hpeak{$pid}{downloc} > ($txstart-$hpeak{$pid}{center})) { # update "downstream of tss" info
		    $hpeak{$pid}{downgene} = $id;
		    $hpeak{$pid}{downloc} = $txstart-$hpeak{$pid}{center};
		}
	    } elsif ($hpeak{$pid}{center} > $txend) {		# peak is upstream of gene
		if ($hpeak{$pid}{uploc} > ($hpeak{$pid}{center}-$txend)) { # update "upstream of gene" info
		    $hpeak{$pid}{upgene} = $id;
		    $hpeak{$pid}{uploc} = $hpeak{$pid}{center}-$txend;
		}
	    }
	}			# end of sense check
    }				# end of for-loop over peaks
}
close(FP);
print STDERR "done. \n";


for my $pid (keys %hpeak) {
    my %hgene=();
    foreach my $gene(@{$hpeak{$pid}{intx}}) {
	$hgene{$id2gene{$gene}}=1;
    }
    my $ingene;
    my $intxlist;
    if (keys %hgene) {
	foreach my $x (keys %hgene) {
	    $ingene = $x;
	    $intxlist = join(",",@{$hpeak{$pid}{intx}});
	}
    } else {
	$ingene = "N/A";
	$intxlist = "N/A";
    }
    print STDOUT join("\t",
		      $pid,
		      $hpeak{$pid}{score},
		      $ingene,
		      $genelen{$ingene},
		      $id2gene{$hpeak{$pid}{upgene}},
		      -$hpeak{$pid}{uploc},
		      $id2gene{$hpeak{$pid}{downgene}},
		      $hpeak{$pid}{downloc},
		      $intxlist
	)."\n";
}
print STDERR "done.\n";


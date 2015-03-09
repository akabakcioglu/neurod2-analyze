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
    $hpeak{$id}{uptssgene} = "N/A";
    $hpeak{$id}{uptssloc} = $VERYLARGENUMBER;
    $hpeak{$id}{downtssgene} = "N/A";
    $hpeak{$id}{downtssloc} = $VERYLARGENUMBER;
    $hpeak{$id}{downgene} = "N/A";
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
	my ($id,$genename) = split/\t/;
	$id2gene{$id}=$genename;
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
    next if ($id2gene{$id} =~ /^Gm|Rik$/); # skip transcripts starting with "Gm.."
    my $glen = $txend-$txstart;
    my $genename = $id2gene{$id};
    $genelen{$genename} = $glen unless (defined $genelen{$genename} && $genelen{$genename} >= $glen);

    for my $pid (keys %{$hchr{$chr}}) { # check if any peak falls here
	my $ctr = $hpeak{$pid}{center};
	if ($ctr > $txstart && $ctr < $txend) { # peak center in transciption region
	    push @{$hpeak{$pid}{intx}},$id;
	}
	# check up/down stream of tx-start
	if ($sense eq "+") { # upstream is smaller than $txstart
	    if ($ctr < $txstart) { # peak is upstream of tss
		if ($hpeak{$pid}{uptssloc} > ($txstart-$hpeak{$pid}{center})) { # update "upstream of tss" info
		    $hpeak{$pid}{uptssgene} = $id;
		    $hpeak{$pid}{uptssloc} = $txstart-$hpeak{$pid}{center};
		}
	    } else {		# peak is downstream of tss
		if ($hpeak{$pid}{downtssloc} > ($hpeak{$pid}{center}-$txstart)) { # update "downstream of tss" info
		    $hpeak{$pid}{downtssgene} = $id;
		    $hpeak{$pid}{downtssloc} = $hpeak{$pid}{center}-$txstart;
		}
		if (($hpeak{$pid}{center}-$txend) > 0) { # peak is downstream of gene-end
		    if ($hpeak{$pid}{downloc} > ($hpeak{$pid}{center}-$txend)) { # update "downstream of gene-end" info
			$hpeak{$pid}{downgene} = $id;
			$hpeak{$pid}{downloc} = $hpeak{$pid}{center}-$txend;
		    }
		}
	    }
	} else {	      # if sense is "-" switch "up" and "down"
	    if ($hpeak{$pid}{center} < $txend) { # peak is downstream of tss
		if ($hpeak{$pid}{downtssloc} > ($txend-$hpeak{$pid}{center})) { # update "downstream of tss" info
		    $hpeak{$pid}{downtssgene} = $id;
		    $hpeak{$pid}{downtssloc} = $txend-$hpeak{$pid}{center};
		}
		
		if ($hpeak{$pid}{center} < $txstart) { # peak is downstream of gene end
		    if ($hpeak{$pid}{downloc} > ($txstart-$hpeak{$pid}{center})) { # update "downstream of gene-end" info
			$hpeak{$pid}{downgene} = $id;
			$hpeak{$pid}{downloc} = $txstart-$hpeak{$pid}{center};
		    }
		}		
	    } else {		# peak is upstream of tss
		if ($hpeak{$pid}{uptssloc} > ($hpeak{$pid}{center}-$txend)) { # update "upstream of gene-end" info
		    $hpeak{$pid}{uptssgene} = $id;
		    $hpeak{$pid}{uptssloc} = $hpeak{$pid}{center}-$txend;
		}
	    }
	}			# end of sense check
    }				# end of for-loop over peaks
}
close(FP);
print STDERR "done.\n";


print STDERR "printing data: merge_id / score / in_gene1,in_gene2,.. / upstream_of / upstream_by / downstream_of(tss) / downstream_by(tss) / ensid1,ensid2,..\n";

for my $pid (keys %hpeak) {
    my %hgene=();
    foreach my $id(@{$hpeak{$pid}{intx}}) {
	$hgene{$id2gene{$id}}=1;
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
		      $id2gene{$hpeak{$pid}{uptssgene}},
		      -$hpeak{$pid}{uptssloc},
		      $id2gene{$hpeak{$pid}{downtssgene}},
		      $hpeak{$pid}{downtssloc},
		      $id2gene{$hpeak{$pid}{downgene}},
		      $hpeak{$pid}{downloc},
		      $intxlist
	)."\n";
}
print STDERR "done.\n";


#/usr/bin/perl -w
use strict;

$|=1;				# flush print commands

my $VERYLARGENUMBER=100000000;
my %hpeak=();
my %hchr=();

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
my $filename = "ensNCrna.txt";

open(FP,$pathname.$filename) or die "Can't open $filename.\n";

print STDERR "Reading ncrna transcript data ($filename) and locating peaks ..\n";
my %hgene = ();
my %id2gene=();
$id2gene{"N/A"}="N/A";
## remove the header if there is one
<FP>; unless (/^\#/) {close(FP); open(FP,$pathname.$filename) or die "Can't open $filename.\n";} # reset pointer

my %rnatype=();
my ($chr,$txstart,$txend,$ncid,$rnatype,$sense,$txid);

while(<FP>) {
    chomp;
    ($chr,$txstart,$txend,$ncid,$rnatype,$sense,$txid) = split/\t/;
    $rnatype{$ncid}=$rnatype;
    for my $pid (keys %{$hchr{$chr}}) { # check if any peak falls here
	my $ctr = $hpeak{$pid}{center};
	if ($ctr > $txstart && $ctr < $txend) { # peak center in transciption region
	    $hpeak{$pid}{ncrna}{$ncid}=1;
	}
	# check up/down stream of tx-start
	if ($sense eq "+") { # upstream is smaller than $txstart
	    if ($ctr < $txstart) { # peak is upstream of tss
		if ($hpeak{$pid}{uploc} > ($txstart-$hpeak{$pid}{center})) { # update "upstream of tss" info
		    $hpeak{$pid}{upgene} = $ncid;
		    $hpeak{$pid}{uploc} = $txstart-$hpeak{$pid}{center};
		    }
	    } else {		# peak is downstream of tss
		if ($hpeak{$pid}{downloc} > ($hpeak{$pid}{center}-$txstart)) { # update "downstream of tss" info
		    $hpeak{$pid}{downgene} = $ncid;
		    $hpeak{$pid}{downloc} = $hpeak{$pid}{center}-$txstart;
		    }
	    }
	} else {	      # if sense is "-" switch "up" and "down"
	    if ($hpeak{$pid}{center} < $txend) { # peak is downstream of tss
		if ($hpeak{$pid}{downloc} > ($txend-$hpeak{$pid}{center})) { # update "downstream of tss" info
		    $hpeak{$pid}{downgene} = $ncid;
		    $hpeak{$pid}{downloc} = $txend-$hpeak{$pid}{center};
		}
	    } else {		# peak is upstream of gene
		if ($hpeak{$pid}{uploc} > ($hpeak{$pid}{center}-$txend)) { # update "upstream of gene" info
		    $hpeak{$pid}{upgene} = $ncid;
		    $hpeak{$pid}{uploc} = $hpeak{$pid}{center}-$txend;
		}
	    }
	}			# end of sense check
    }				# end of for-loop over peaks
}
close(FP);
print STDERR "done. Converting ensembl IDs to gene names..\n";

my $filename = "/home/alkan/c/gulayse/neurod2/mouse_data/ucsc/ensemblToGeneName.txt";
open(FP,$filename) or die "Can't open ensemblToGeneName.txt";
while(<FP>) {
    chomp;
    my ($id,$geneid) = split/\t/;
    $id2gene{$id}=$geneid;
}

print STDERR "printing data: merge_id / score / in_gene1,in_gene2,.. / upstream_of / upstream_by / downstream_of(tss) / downstream_by(tss) / ensid1,ensid2,..\n";

for my $pid (keys %hpeak) {
    # my %hgene=();
    # foreach my $gene(@{$hpeak{$pid}{intx}}) {
    # 	$hgene{$id2gene{$gene}}=1;
    # }
    my $ingenelist;
    my $intxlist;
    # if (keys %hgene) {
    # 	$ingenelist = join(",",sort keys %hgene);
    # 	$intxlist = join(",",@{$hpeak{$pid}{intx}});
    if (defined $hpeak{$pid}{ncrna}) {
	my $typ;
	foreach my $nc (keys %{$hpeak{$pid}{ncrna}}) {
	    $ingenelist .= $nc.",";
	    $typ = $rnatype{$nc};
	}
	$ingenelist .= "$typ";
    } else {
	$ingenelist = "N/A";
	$intxlist = "N/A";
    }
    print STDOUT join("\t",
		      $pid,
		      $hpeak{$pid}{score},
		      $ingenelist,
		      $hpeak{$pid}{upgene},
		      -$hpeak{$pid}{uploc},
		      $hpeak{$pid}{downgene},
		      $hpeak{$pid}{downloc}
	)."\n";
}
print STDERR "done.\n";


#/usr/bin/perl -w
use strict;
use List::MoreUtils qw(uniq);

$|=1;				# flush print commands

my $VERYLARGENUMBER=100000000;
my %hpeak=();
my %hchr=();

my $exppath = shift @ARGV;
my $filetype = shift @ARGV;
die "Run as: > perl $0.pl (path_to_exp_data) (ens|refseq) < inputfile > outfile. Specified type must be either ens or refseq!\n" 
    unless ($filetype =~ /ens|refseq/);


print STDERR "Reading merge.txt .. ";
my $filename = "merge.txt";
open(FP,$exppath."/".$filename) or die "Can't open $filename";
while (<FP>){			# read one sequence at a time
    chomp;
    my ($pid,$nmerge,$score,$chr,$start,$end,$lmerge,$mstart,$mend,@mergelist) = split/\t/;
    $hpeak{$pid}{score} = $score;
    $hpeak{$pid}{chr} = $chr;
    $hpeak{$pid}{start} = $start;
    $hpeak{$pid}{end} = $end;
    $hpeak{$pid}{center} = int(0.5*($start+$end));
    $hpeak{$pid}{region} = "N/A";
    $hchr{$chr}{$pid}=1;
}
close(FP);
my $npeak = keys %hpeak;
print STDERR "finished reading $npeak peaks.\n";




if ($filetype =~ /ens/) {
    $filename = "context_ens.txt";
} elsif ($filetype =~ /refseq/) {
    $filename = "context_refseq.txt";
}

print STDERR "Reading $filename (to get peaks located in genes).\n";
open(FP,$exppath."/".$filename) or die "Can't open $filename";

my %hpeaksintx=();

while (<FP>){			# read one sequence at a time
    chomp;
    my ($pid,$score,$genes,$upgene,$upby,$downgene,$downby,$ids) = split/\t/;
    my @genelist = split/,/,$genes;
    my @idlist = split/,/,$ids;
    @{$hpeak{$pid}{genenames}} = @genelist;
    @{$hpeak{$pid}{geneids}} = @idlist;
    foreach my $tx (@idlist) {
	push @{$hpeaksintx{$tx}},$pid;
    }
}
close(FP);
print STDERR "finished reading $filename. Found " . scalar(keys %hpeaksintx) . " transcripts occupied by peaks\n";



my $mousepath = "/home/alkan/c/gulayse/neurod2/mouse_data/ucsc/";

my %htx=();

# Read transcript starts/ends
if ($filetype =~ /ens/) {
    $filename = "ensGene_uniq.txt";
    print STDERR "Reading transcript data ($filename).\n";
    open(FP,$mousepath.$filename) or die "Can't open $filename";
    while(<FP>) {
	chomp;
	my ($chr,$start,$end,$txid,@rest) = split/\t/;
	$htx{$txid}{start} = $start;
	$htx{$txid}{end} = $end;
	$htx{$txid}{chr} = $chr;
    }
} elsif ($filetype =~ /refseq/) {
    $filename = "RefSeq_refgene_uniq.txt";
    print STDERR "Reading transcript data ($filename).\n";
    open(FP,$mousepath.$filename) or die "Can't open $filename";
    <FP>;			# header
    while(<FP>) {
	chomp;
	my ($bin,$txid,$chr,$sense,$start,$end,@rest) = split/\t/;
## CAREFUL, RefSeq geneid is not unique!! May map to multiple positions/chromosomes..
	unless (defined $htx{$txid} && # dont record position if not a new id and
		($end-$start) < ($htx{$txid}{end}-$htx{$txid}{start})) { # the new location is shorter..
		$htx{$txid}{start} = $start;
		$htx{$txid}{end} = $end;
		$htx{$txid}{chr} = $chr;
	} # this way, the probability ($p below) will be well behaved.
    }
}

close(FP);
print STDERR "Finished reading transcript data ($filename).\n";




# Read each chromosome in mm10 and locate the peaks within genes
my $mousepath = "/home/alkan/c/gulayse/neurod2/mouse_data/ucsc/";

if ($filetype =~ /ens/) {
    $filename = "ensRegion_uniq.txt";
} elsif ($filetype =~ /refseq/) {
    $filename = "RefSeq_Region_uniq.txt";
}

open(FP,$mousepath.$filename) or die "Can't open $filename";

print STDERR "Reading intron/exon location data ($filename) and locating peaks ..\n";
my %hreglen=();
my %hreglist=();
my %hregcount=();
my %hregscore=();

my %hmean=(); # calculates mean hit counts based on region lengths
my %hvar=(); # calculates variance of hit counts

while(<FP>) {
    chomp;
    my ($chr,$start,$end,$txid,$sense,$regtype,$order,$regionid) = split/\t/;

    $hreglen{$regtype}{$order} += $end-$start; 
    $hregcount{$regtype}{$order} += 1; # initially zero by default
    @{$hreglist{$regtype}{$order}} = () unless (defined $hreglist{$regtype}{$order});

    for my $pid (@{$hpeaksintx{$txid}}) { # check if any peak falls here
	my $ctr = $hpeak{$pid}{center};
	my $p = (1.0*($end-$start))/($htx{$txid}{end}-$htx{$txid}{start});
	$hmean{$regtype}{$order} += $p; # mean of a coin flip with prob p
	$hvar{$regtype}{$order} += $p*(1-$p); # variance of a coin flip with prob p
	if ($ctr >= $start && $ctr <= $end) { # peak in region
	    push @{$hreglist{$regtype}{$order}},$pid;
	    $hregscore{$regtype}{$order} += $hpeak{$pid}{score};
	}
    }				# end of for-loop over peaks
}
close(FP);


## print data
for my $reg (keys %hreglist) {
    for my $ord (keys %{$hreglist{$reg}}) {
	my $npeaks = scalar @{$hreglist{$reg}{$ord}};
	
	if (			#$hmean{$reg}{$ord}>1  && 
	    $npeaks > 0) {
	    my %hreggenes=();
	    foreach my $p (@{$hreglist{$reg}{$ord}}) {
		foreach my $g (@{$hpeak{$p}{genenames}}) {
		    $hreggenes{$g} += 1;
		}
	    }
	    print join("\t",
		       $reg,
		       $ord,
		       $npeaks,
		       $hreglen{$reg}{$ord}/$hregcount{$reg}{$ord},
		       ($npeaks-$hmean{$reg}{$ord})/sqrt($hvar{$reg}{$ord}), # my z-score 
		       $hregscore{$reg}{$ord}/$npeaks,
#		       join(",",@{$hreglist{$reg}{$ord}})
		       join(",",keys %hreggenes),
		       join(",",uniq @{$hreglist{$reg}{$ord}})
		)."\n";
	} elsif ($hvar{$reg}{$ord} < 0) {
	    print STDERR join("\t",$hmean{$reg}{$ord},$hvar{$reg}{$ord},$npeaks)."***\n";
	}
    }
}

print STDERR "done.\n";

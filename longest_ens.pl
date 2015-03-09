#/usr/bin/perl -w
use strict;

$|=1;				# flush print commands

my $VERYLARGENUMBER=100000000;
my %hchr=();

my $filetype = shift @ARGV;
die "Run as: > perl $0.pl (ens|refseq) < inputfile > outfile. Specified type must be either ens or refseq!\n" 
    unless ($filetype =~ /ens|refseq/);

# Read each chromosome in mm10 and locate the peaks within genes
my $pathname = "/home/alkan/c/gulayse/neurod2/mouse_data/ucsc/";
my $filename;

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
close(FP);

print STDERR "Reading transcript data ($filename) and locating peaks ..\n";
if ($filetype =~ /ens/) {
    $filename = "ensGene.txt";
} elsif ($filetype =~ /refseq/) {
    $filename = "RefSeq_refgene.txt";
}
open(FP,$pathname.$filename) or die "Can't open $filename.\n";

## remove the header if there is one
<FP>; unless (/^\#/) {close(FP); open(FP,$pathname.$filename) or die "Can't open $filename.\n";} # reset pointer

my @unused=();
my ($chr,$txstart,$txend,$id,$sense,$bin,$cdsstart,$cdsend,$exoncnt,$exonstarts,$exonends,$b,$genename);

my %hlongest_tx_ofgene=();
my %hgenelength=();

while(<FP>) {
    chomp;
    if ($filetype =~ /ens/) {	# checking again and again is stupid, but let's not worry about it
	($chr,$txstart,$txend,$id,$a,$sense,@unused) = split/\t/;
    } elsif ($filetype =~ /refseq/) {
	($bin,$id,$chr,$sense,$txstart,$txend,$cdsstart,$cdsend,$exoncnt,$exonstarts,$exonends,$b,$genename,@unused) = split/\t/;
	$id2gene{$id}=$genename;
    }
    my $txlen = $txend-$txstart;
    unless (defined $hgenelength{$id2gene{$id}} 
	    and $hgenelength{$id2gene{$id}} > $txlen) {
	$hlongest_tx_ofgene{$id2gene{$id}} = $id;
	$hgenelength{$id2gene{$id}} = $txlen;
    }
}

my %hlongest_tx=();
foreach my $id (keys %hlongest_tx_ofgene) {
    $hlongest_tx{$hlongest_tx_ofgene{$id}} = $id;
}

# reopen the file to pick relevant lines.
close(FP); open(FP,$pathname.$filename) or die "Can't open $filename.\n";
## remove the header if there is one
<FP>; unless (/^\#/) {close(FP); open(FP,$pathname.$filename) or die "Can't open $filename.\n";} # reset pointer
while(<FP>) {
    if ($filetype =~ /ens/) {	# checking again and again is stupid, but let's not worry about it
	($chr,$txstart,$txend,$id,$a,$sense,@unused) = split/\t/;
    } elsif ($filetype =~ /refseq/) {
	($bin,$id,$chr,$sense,$txstart,$txend,$cdsstart,$cdsend,$exoncnt,$exonstarts,$exonends,$b,$genename,@unused) = split/\t/;
    }
    print if (defined $hlongest_tx{$id});
}

close(FP);

## also filter these ens tx's from ensRegion.txt

if ($filetype =~ /ens/) {
    $filename = "ensRegion.txt";
    my $outfile = "ensRegion_uniq.txt";
    open(FP,$pathname.$filename) or die "Can't open $filename.\n";
    open(FPOUT,">".$pathname.$outfile) or die "Can't open $filename.\n";
## remove the header if there is one
    while(<FP>) {
	($chr,$txstart,$txend,$id,@unused) = split/\t/;
	print FPOUT if (defined $hlongest_tx{$id});
    }
    close(FP);
}

if ($filetype =~ /refseq/) {
    $filename = "RefSeq_Region.txt";
    my $outfile = "RefSeq_Region_uniq.txt";
    open(FP,$pathname.$filename) or die "Can't open $filename.\n";
    open(FPOUT,">".$pathname.$outfile) or die "Can't open $filename.\n";
## remove the header if there is one
    while(<FP>) {
	($chr,$txstart,$txend,$id,@unused) = split/\t/;
	print FPOUT if (defined $hlongest_tx{$id});
    }
    close(FP);
}

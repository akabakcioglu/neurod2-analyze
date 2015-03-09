#/usr/bin/perl -w
use strict;

$|=1;				# flush print commands

my $filetype = shift @ARGV;	# ens or refseq
die "Run as $0 <ens|refseq> < infile > outfile\n" unless ($filetype =~ /ens|refseq/);

# # Read each chromosome in mm10 and locate the peaks within genes
# my $fname = "ensRegion_clean.txt";
# my $pathname = "/home/alkan/c/gulayse/neurod2/mouse_data/ucsc/";
# open(FP,$pathname.$fname) or die "Can't open $fname";

my %h=();
my %hmaxorder=();

#while(<FP>) {
while(<>) {
    chomp;
    my ($chr,$start,$end,$regionid,$a,$sense) = split/\t/;
    my ($head,$id,$regname,$order,,@rest);
    if ($filetype =~ /ens/) {
	($id,$regname,$order,@rest) = split/_/,$regionid;
    } elsif ($filetype =~ /refseq/) {
	($head,$id,$regname,$order,@rest) = split/_/,$regionid;
	$id = "$head"."_"."$id";
    }

    $h{$regionid}{ensid}=$id;
    $h{$regionid}{chr}=$chr;
    $h{$regionid}{start}=$start;
    $h{$regionid}{end}=$end;
    $h{$regionid}{sense}=$sense;
    $h{$regionid}{reg}=$regname;
    $h{$regionid}{order}=$order;
    $hmaxorder{$id}=$order
	unless (defined $hmaxorder{$id} && $hmaxorder{$id} >= $order);
    if ($h{$regionid}{reg} eq "intron") {
	$h{$regionid}{isexon}=0;
    } else {
	$h{$regionid}{isexon}=1;
    }
}

for my $x (sort keys %h) {
    if ($h{$x}{sense} eq "-") {
	$h{$x}{order} = $hmaxorder{$h{$x}{ensid}} - $h{$x}{order};
	$h{$x}{order} -= 1 unless ($h{$x}{isexon});
    }
}

## correction for utr3 (start counting from the end

for my $x (sort keys %h) {
    if ($h{$x}{reg} eq "utr3") {
	$h{$x}{order} = $hmaxorder{$h{$x}{ensid}} - $h{$x}{order};
    }
}

for my $x (sort keys %h) {
    print join("\t",$h{$x}{chr},$h{$x}{start},$h{$x}{end},$h{$x}{ensid},$h{$x}{sense},$h{$x}{reg},$h{$x}{order},$x)."\n";
}
    
    

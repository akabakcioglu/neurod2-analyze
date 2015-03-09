#/usr/bin/perl -w
use strict;

$|=1;				# flush print commands

my $datapath = "/home/alkan/c/gulayse/neurod2/gulayse_data/E14_merge_new/";
my $filename = "merge.txt";
open(FP,$datapath.$filename) or die "Can't open $filename";
my %hmergectr=();
while (<FP>){			# read one sequence at a time
    chomp;
    my ($id,$nmerge,$score,$chr,$start,$end,$lmerge,$mstart,$mend,@mergelist) = split/\t/;
    $hmergectr{$id} = int(0.5*($start+$end));
}
close(FP);
print STDERR "finished reading $filename.\n";

my $mousepath = "/home/alkan/c/gulayse/neurod2/mouse_data/ucsc/";
my $filename = "ensRegion.txt";
open(FP,$mousepath.$filename) or die "Can't open $filename";
my %hstart=();
my %hend=();
my %htype=();
while(<FP>) {
    chomp;
    my ($chr,$start,$end,$txid,$sense,$regtype,$order,$regionid) = split/\t/;
    $hstart{$txid}{$regionid}=$start;
    $hend{$txid}{$regionid}=$end;
    $htype{$txid}{$regionid}=$regtype;
}
close(FP);
print STDERR "finished reading $filename.\n";


my %hscore=();
my %hlen=();

my $margin = 10000;
my $fmargin = 10000.0;

my $cnt = 0;
while (<>){			# read context_ens.txt
    $cnt += 1;
    print STDERR "." if ($cnt % 100 == 0);
    chomp;
    my ($id,$score,$ingene,$inlen,$uptssgene,$uptsslen,$downtssgene,$downtsslen,$downgene,$downlen,$ensliststring) = split/\t/;
    my @enslist = split/,/,$ensliststring;
    my %h=();
    if ($ingene !~ /N\/A/ && $inlen =~ /[0-9]*/) {
	$hlen{$ingene}=$inlen;
	my $sc = $fmargin/$inlen;
	my $inintron = 0;
	foreach my $tx (@enslist) { # consider each tx in which this peak appears
	    my @reglist = sort {$hstart{$tx}{$a}<=>$hstart{$tx}{$b}} keys %{$hstart{$tx}}; # sort regions of tx according to start pos
	    my $reg;
	    do {
		$reg = shift @reglist;
	    } while ($hmergectr{$id} > $hend{$tx}{$reg});
	    die "end is after ctr but start is not before: $id\t$tx\t$hmergectr{$id}\t$reg\t$hstart{$tx}{$reg}\t$hend{$tx}{$reg}!!\n" unless ($hmergectr{$id} >= $hstart{$tx}{$reg});
	    $inintron = 1 if ($htype{$tx}{$reg} =~ /intron/);
	}	    
	$h{$ingene} = $sc>1.0?1:$sc if ($inintron == 1);
#	$h{$ingene} = $sc 

    }
    if (-$uptsslen <= $margin) {	# $uplen is negative
	$h{$uptssgene} = 1 unless (defined $h{$uptssgene} && $h{$uptssgene} > 1);
    }
    if ($downlen <= $margin) {	# $uplen is negative
	$h{$downgene} = 1 unless (defined $h{$downgene} && $h{$downgene} > 1);
    }
    for my $x (keys %h) {
	$hscore{$x} += $h{$x};
    }
}    

for my $x (reverse sort {$hscore{$a}<=>$hscore{$b}} keys %hscore) {
    print "$x\t$hscore{$x}\t$hlen{$x}\n";
}



#/usr/bin/perl -w
use strict;
use List::Util qw[min max];

my $VERYLARGE = 1e15;
my $score;
my %hpeaks=();
my %hchr=();
my %hstart=();
my %hend=();
my %hscore=();
my $cnt = 0;

my ($datapath, $overlap_thr, @filelist) = @ARGV;


## READ PEAK DATA
foreach my $fname (@filelist) {
    my $filename = $datapath . $fname;
    open(FP,$filename) or die "Can't open".$filename;
    my $id_cnt=0;
    while (<FP>) {	       # 4-column input file: chr start end score
	chomp;
	my ($chr,$start,$end,$length,$summit,$ntags,$score,$fold,$fdr) = split/\t/;
#	my ($chr,$start,$end,$peakid,$score) = split/\t/;
	$id_cnt += 1;
	my $id = $fname."_tag_".$id_cnt;
	$hstart{$id}=$start;
	$hend{$id}=$end;
	$hchr{$id}=$chr;
	$hscore{$id}=$score;
	$hpeaks{$chr}{$id}=$score;
    }
    close(FP);
}

## OVERLAP CALCULATOR BETWEEN TWO PEAKS
sub overlap {
    my ($p1,$p2) = @_;
    my $overlap = 0;
    my ($s1,$e1,$s2,$e2) = ($hstart{$p1},$hend{$p1},$hstart{$p2},$hend{$p2});
    return $overlap if ($hchr{$p1} != $hchr{$p2});
    return $overlap if ($e1 < $s2 or $e2 < $s1);
    my $l1 = $e1 - $s1;
    my $l2 = $e2 - $s2;
    my $lmin = $l1 < $l2 ? $l1 : $l2;
    my $st = $s2 > $s1 ? $s2 : $s1; # $st is the larger start position
    my $end = $e2 > $e1 ? $e1 : $e2; # $end is the smaller end position
    $overlap = $end - $st;
    return $overlap/$lmin; # we want 100% overlap if one peak is inside the other, hence lmin..
}

print STDERR "Finished reading peak positions.\n";

## START CALCULATING MERGES WITH GIVEN OVERLAP THRESHOLD. THAT MUCH OVERLAP
## MUST EXIST BETWEEN EACH PAIR OF THE MERGED PEAK SET.
my %hmerge = ();
my %hmergeid = ();
my %hmchr = ();
my %hchrm = ();
my $nextid=1;


## RECORD MERGE START/END/SCORE/INTERSECTION_REGION INFO

my %hmstart=();			# merge starts
my %hmend=();			# merge ends
my %hmscore=();			# merge scores

##

sub initmerge {
# init start/end/score
    my ($m,$p,$c) = @_;		# read merge and peak and $chr ids
    $hmerge{$m}{$p}=1;
    $hmergeid{$p} = $m;
    ## update start/end/score
    $hmscore{$m} = $hscore{$p};
    $hmstart{$m} = $hstart{$p};
    $hmend{$m} = $hend{$p};
    # 
    $hchrm{$c}{$m}=1;
    $hmchr{$m}=$c;
}

sub removemerge {
    my $m = shift;
    delete $hchrm{$hmchr{$m}} if defined $hchrm{$hmchr{$m}};
    delete $hmchr{$m} if defined $hmchr{$m};
    delete $hmend{$m} if defined $hmend{$m};
    delete $hmstart{$m} if defined $hmstart{$m};
    delete $hmscore{$m} if defined $hmscore{$m};
    delete $hmerge{$m} if defined $hmerge{$m};
}

sub addpeaktomerge {
# update start/end/score
    my ($m,$p) = @_;		# read merge and peak ids
    unless (defined $hmerge{$m}{$p}) {
	$hmerge{$m}{$p}=1;
	$hmergeid{$p} = $m;
	$hmstart{$m} = $hstart{$p} if ($hstart{$p}<$hmstart{$m});
	$hmend{$m} = $hend{$p} if ($hend{$p}>$hmend{$m});
	$hmscore{$m} += $hscore{$p};
    }
}

sub getoverlap {
    # returns the min,max overlap with the given peak and 
    # max overlapping peak in the given merge
    my ($m,$p) = @_;
    my $min = 1.0;
    my $max = 0.0;
    my $max_id;
    foreach my $pid (keys %{$hmerge{$m}}) {
	my $o = overlap($pid,$p);
	$min = $o if ($o < $min);
	if ($o > $max) {
	    $max = $o;
	    $max_id = $pid;
	}
    }
    return ($min,$max,$max_id);
}


sub mergemerges {
    # merge two merges into the first one
    my ($mkeep,$mremove) = @_;
    for my $p (keys %{$hmerge{$mremove}}) {
	addpeaktomerge($mkeep,$p) unless (defined $hmerge{$mkeep}{$p});
    }
    removemerge($mremove);
}


print STDERR "Merging peaks..\n";
foreach my $chr (sort keys %hpeaks) {
    next if ($chr =~ /random|chrUn/);
    print STDERR " $chr ";
    my @idlist = sort {$hstart{$a}<=>$hstart{$b}} keys %{$hpeaks{$chr}};
    print STDERR "includes ".scalar(keys %{$hpeaks{$chr}})." peaks\n";

    foreach my $id (@idlist) { # sort wrt start pos in increasing order

	my @mids = sort {$hmend{$a}<=>$hmend{$b}} keys %{$hchrm{$chr}}; 

	my $mergeid = pop @mids; # start from the last merge-end pos and proceed backwards

	while ($mergeid and $hmend{$mergeid} >= $hstart{$id}) { # run over merges, sorted
	    die "ERROR!!! in merge.pl\n" if ($hmstart{$mergeid} > $hstart{$id}); # with this sorting, shouldn't happen
	    my ($minoverlap,$maxoverlap,$maxoverlap_id) = getoverlap($mergeid,$id);
	    if ($maxoverlap >= $overlap_thr) { # add id to merge, good overlap with at least one peak
		unless (defined ($hmergeid{$id})) { # first time an overlapping merge is found
		    addpeaktomerge($mergeid,$id);
		} else {       # already merged $id. now merge merges!
		    mergemerges($mergeid,$hmergeid{$id}); # order matters. keep the new one, remove old merge.
		}
	    }
	    $mergeid = pop @mids;
	}		       # further merges do not overlap with id
	unless (defined $hmergeid{$id}) {
	    initmerge($nextid,$id,$chr);
	    $nextid += 1;
	}
    }
}


foreach my $mergeid (reverse sort{$hmscore{$a}<=>$hmscore{$b}} keys %hmerge) { # list peaks in each merge
    if (scalar (keys %{$hmerge{$mergeid}}) > 1) {
	print join("\t",$mergeid,scalar (keys %{$hmerge{$mergeid}}),$hmscore{$mergeid},$hmchr{$mergeid},$hmstart{$mergeid},$hmend{$mergeid},$hmend{$mergeid}-$hmstart{$mergeid},$hmstart{$mergeid},$hmend{$mergeid},keys %{$hmerge{$mergeid}})."\n";
    }
}

print STDERR "Done.\n";

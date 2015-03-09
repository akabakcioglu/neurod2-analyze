#/usr/bin/perl -w
use strict;

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

print STDERR "Merging peaks..\n";
foreach my $chr (sort keys %hpeaks) {
    next if ($chr =~ /random|chrUn/);
    print STDERR " $chr ";
    my @idlist = sort {$hstart{$a}<=>$hstart{$b}} keys %{$hpeaks{$chr}};
    print STDERR "includes ".scalar(keys %{$hpeaks{$chr}})." peaks\n";
    foreach my $id (@idlist) { # sort wrt start pos in increasing order
	my @ids = @idlist;     # temporary list of ids to scan
	my $newid;
	do {
	    $newid = shift @ids;	# pick first

	    if ($hstart{$newid} < $hend{$id} and
		$hstart{$id} < $hend{$newid}) { # there's overlap

		unless ($newid eq $id) {
		
		    foreach my $mergeid (keys %{$hmergeid{$newid}}) {
			my $minoverlap = 1.0;
			foreach my $pid (keys %{$hmerge{$mergeid}}) {
			    my $newoverlap = overlap($pid,$id);
			    $minoverlap = $newoverlap if ($newoverlap < $minoverlap);
			}
			if ($minoverlap >= $overlap_thr) { # add id to merge
			    $hmerge{$mergeid}{$id} = 1;
			    $hmergeid{$id}{$mergeid} = 1;
			}
		    }
		}
	    }
	} while (@ids and $hstart{$newid} <= $hend{$id}); # further newids do not overlap with id
	unless (defined $hmergeid{$id}) {
	    $hmerge{$nextid}{$id}=1;
	    $hmergeid{$id}{$nextid}=1;
	    $hchrm{$chr}{$nextid}=1;
	    $hmchr{$nextid}=$chr;
	    $nextid += 1;
	}
    }
}

print STDERR "\nFinished merging peaks. Now extracting intersecting segments..\n";

## RECORD MERGE START/END/SCORE/INTERSECTION_REGION INFO

my %hmstart=();			# merge starts
my %hmend=();			# merge ends
my %hmscore=();			# merge scores
my %hmintst=();			# merge intersection start
my %hmintend=();		# merge intersection end

foreach my $mergeid (sort{$a<=>$b} keys %hmerge) { # list peaks in each merge
    $hmstart{$mergeid}=$VERYLARGE;
    $hmend{$mergeid}=0;
    $hmscore{$mergeid}=0;
    $hmintst{$mergeid}=0;	# intersection region of merge - start
    $hmintend{$mergeid}=$VERYLARGE; # intersection region of merge - end
    foreach my $id (keys %{$hmerge{$mergeid}}) {
	$hmstart{$mergeid}=$hstart{$id} if ($hmstart{$mergeid} > $hstart{$id});
	$hmend{$mergeid}=$hend{$id} if ($hmend{$mergeid} < $hend{$id});
	$hmscore{$mergeid} += $hscore{$id};
	$hmintst{$mergeid}=$hstart{$id} if ($hmintst{$mergeid} < $hstart{$id});
	$hmintend{$mergeid}=$hend{$id} if ($hmintend{$mergeid} > $hend{$id});
    }
}


foreach my $mergeid (reverse sort{$hmscore{$a}<=>$hmscore{$b}} keys %hmerge) { # list peaks in each merge
    if (scalar (keys %{$hmerge{$mergeid}}) > 1) {
	print join("\t",$mergeid,scalar (keys %{$hmerge{$mergeid}}),$hmscore{$mergeid},$hmchr{$mergeid},$hmintst{$mergeid},$hmintend{$mergeid},$hmintend{$mergeid}-$hmintst{$mergeid}+1,$hmstart{$mergeid},$hmend{$mergeid},keys %{$hmerge{$mergeid}})."\n";
    }
}

print STDERR "Done.\n";

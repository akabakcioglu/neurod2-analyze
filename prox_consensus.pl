#/usr/bin/perl -w
use strict;
my $id;
my $dir = shift @ARGV;
#my $cs = "ca[gt][ca]tg[gt]";
my $cs = "cagatgg";
my $cs_rev = "ccatctg";
while (<>) {
    chomp;
    if (/^>(.*)$/) {
	$id = $1;
    } elsif ($dir =~ /\-/) { # get upstream of the cs
	while (/(.{0,100})$cs/ig) {
	    my $prevseq = $1;
	    print STDOUT ">$id\n$prevseq\n" if ($prevseq);
	}
	while (/$cs_rev(.{0,100})/ig) {
	    my $prevseq = $1;
	    print STDOUT ">$id"."_rev\n$prevseq\n" if ($prevseq);
	}
    } elsif ($dir =~ /\+/) {	# get downstream of the cs
	while (/$cs(.{0,100})/ig) {
	    my $postseq = $1;
	    print STDOUT ">$id\n$postseq\n" if ($postseq);
	}
	while (/(.{0,100})$cs_rev/ig) {
	    my $postseq = $1;
	    print STDOUT ">$id"."_rev\n$postseq\n" if ($postseq);
	}
    } elsif ($dir =~ /0/) {	# get both sides of the cs
	while (/(.{100}$cs.{100})/ig) {
	    my $centerseq = $1;
	    print STDOUT ">$id\n$centerseq\n" if ($centerseq);
	}
    }
}

	

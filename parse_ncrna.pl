#/usr/bin/perl -w
use strict;

$|=1;				# flush print commands

my $infilename = "/home/alkan/c/gulayse/neurod2/mouse_data/ucsc/Mus_musculus.GRCm38.75.ncrna.fa";
my $logfilename = "/home/alkan/c/gulayse/neurod2/analyze/parse_ncrna.log";
open(FP,$infilename) or die "cannot open $infilename\n";
open(FPLOG,">".$logfilename) or die "cannot open $logfilename\n";
while(<FP>) {
    my $line = $_;
    chomp($line);
    if($line =~ /^>(.*)$/) {		# new gene
	my ($txid,$ncrnainfo,$locinfo,$geneid,$genebiotype,$txbiotype) = split/ /,$1;
	my ($txt_ncrna,$ncrna) = split/:/,$ncrnainfo;
	my ($txt_chr,$txt_GRCm38,$chr,$start,$end,$sense) = split/:/,$locinfo;
	my ($txt_gene,$gene) = split/:/,$geneid;
	my ($txt_gene_biotype,$gene_biotype) = split/:/,$genebiotype;
	my ($txt_tx_biotype,$tx_biotype) = split/:/,$txbiotype;

	unless ($txid =~ /^ENSMUS/ &&
		$txt_chr eq "chromosome" &&
		$txt_GRCm38 eq "GRCm38" &&
		$chr =~ /^[0-9MXY]*$/ &&
		$start =~ /^[0-9]*$/ &&
		$end =~ /^[0-9]*$/ &&
		($sense == 1 || $sense == -1) &&
		$txt_gene eq "gene" &&
		$gene =~ /^ENSMUS/ &&
		$txt_gene_biotype eq "gene_biotype" && ($ncrna eq $gene_biotype) &&
		$txt_tx_biotype eq "transcript_biotype" && ($ncrna eq $tx_biotype)) {
	    print FPLOG "*** >$line\n";
	    print FPLOG "txid error\n" unless ($txid =~ /^ENSMUS/);
	    print FPLOG "txt_chr error\n" unless ($txt_chr eq "chromosome");
	    print FPLOG "txt_GRCm38 error\n" unless ($txt_GRCm38 eq "GRCm38");
	    print FPLOG "chr error\n" unless ($chr =~ /^[0-9MXY]*$/);
	    print FPLOG "start error\n" unless ($start =~ /^[0-9]*$/);
	    print FPLOG "end error\n" unless ($end =~ /^[0-9]*$/);
	    print FPLOG "sense error\n" unless (($sense == 1 || $sense == -1));
	    print FPLOG "txt_gene error\n" unless ($txt_gene eq "gene");
	    print FPLOG "gene error\n" unless ($gene =~ /^ENSMUS/);
	    print FPLOG "txt_gene_biotype error\n" unless ($txt_gene_biotype eq "gene_biotype" && ($ncrna eq $gene_biotype));
	    print FPLOG "txt_tx_biotype error\n" unless ($txt_tx_biotype eq "transcript_biotype" && ($ncrna eq $tx_biotype));
	} else {	# a legit line
	    my $sign;
	    $sign = "+" if ($sense == 1);
	    $sign = "-" if ($sense == -1);
	    print STDOUT join("\t","chr".$chr,$start,$end,$gene,$ncrna,$sign,$txid),"\n";
	}
    }
}    
close(FP);
close(FPLOG);	    
    
    

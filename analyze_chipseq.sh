## THIS VERSION FOLLOWS THE REFERENCE:
## http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1003342#abstract0
## ----------------------------------------------------------------------------------------

## all steps of the analysis are combined here. The input is a list of
## peak files (generated with GFP control)
DATA_PATH="../gulayse_data/P0_plos"
## find regions where peaks significantly overlap
perl merge.pl $DATA_PATH/ 0.50 R1.txt R2.txt M.txt > $DATA_PATH/merge.txt
perl -e 'while(<>){chomp;@a=split/\t/;print join("\t",$a[3],$a[4],$a[5],$a[0],$a[2])."\n";}' < $DATA_PATH/merge.txt > $DATA_PATH/merge.bed
## get the sequences of the overlaps using BioPerl
perl merge_seq.pl < $DATA_PATH/merge.txt > $DATA_PATH/merge_seq.txt

## locate overlap regions (tss-relative position, nearby genes) using ensembl data
perl context.pl ens < $DATA_PATH/merge.txt | sort -k2 -n -r > $DATA_PATH/context_ens.txt
## same as above, using refseq data
perl context.pl refseq < $DATA_PATH/merge.txt | sort -k2 -n -r > $DATA_PATH/context_refseq.txt


## PICK ONLY THE LONGEST ENS TRANSCRIPT FOR EACH GENE AND USE THIS WHILE IMPLEMENTING THE PLOS METHOD
#perl longest_ens.pl ens < $DATA_PATH/../../mouse_data/ucsc/ensGene.txt > $DATA_PATH/../../mouse_data/ucsc/ensGene_uniq.txt
#perl longest_ens.pl refseq < $DATA_PATH/../../mouse_data/ucsc/RefSeq_refgene.txt > $DATA_PATH/../../mouse_data/ucsc/RefSeq_refgene_uniq.txt

## Implement the plos-algorithmm using the files above.

## get the global up/down-tss distr and record in plos_tss_distr.txt

perl plos_tss_distr.pl ens < $DATA_PATH/merge.txt > $DATA_PATH/plos_tss_distr_ens.txt
perl plos_tss_distr.pl refseq < $DATA_PATH/merge.txt > $DATA_PATH/plos_tss_distr_refseq.txt

## get plos scores (and distributions along the way)

perl plos_score.pl ens < $DATA_PATH/merge.txt > $DATA_PATH/plos_score_ens.txt
perl plos_score.pl refseq < $DATA_PATH/merge.txt > $DATA_PATH/plos_score_refseq.txt

## binding region analysis

## get binding region preferences (5utr/cds/3utr/intron..) on ensembl data
perl bindingregion_plos_v2.pl $DATA_PATH/ ens | sort -k5 -n -r > $DATA_PATH/bindingregion_plos_ens.txt
perl tag_count_in_regions.pl $DATA_PATH/bindingregion_plos_ens.txt > $DATA_PATH/bindingregion_plos_ens.tagcounts
## same as above, using refseq data
perl bindingregion_plos_v2.pl $DATA_PATH/ refseq | sort -k5 -n -r > $DATA_PATH/bindingregion_plos_refseq.txt
perl tag_count_in_regions.pl $DATA_PATH/bindingregion_plos_refseq.txt > $DATA_PATH/bindingregion_plos_refseq.tagcounts
## remember: utr3_0 refers to the first utr3 in reading direction (as for utr5/cds/intron)

perl -e 'while(<>){chomp;@a=split/\t/;if($a[0]=~/utr5/){print $a[6]."\n";}else{print "0\n";}}' < $DATA_PATH/bindingregion_plos_ens.txt > bindingregion_plos_ens.score.utr5
perl -e 'while(<>){chomp;@a=split/\t/;if($a[0]=~/utr3/){print $a[6]."\n";}else{print "0\n";}}' < $DATA_PATH/bindingregion_plos_ens.txt > bindingregion_plos_ens.score.utr3
perl -e 'while(<>){chomp;@a=split/\t/;if($a[0]=~/cds/){print $a[6]."\n";}else{print "0\n";}}' < $DATA_PATH/bindingregion_plos_ens.txt > bindingregion_plos_ens.score.cds
perl -e 'while(<>){chomp;@a=split/\t/;if($a[0]=~/intron/){print $a[6]."\n";}else{print "0\n";}}' < $DATA_PATH/bindingregion_plos_ens.txt > bindingregion_plos_ens.score.intron

use Bio::Perl;
use Bio::DB::Fasta;

# READS A FILE WITH A SEQUENCE PER LINE and CHECK IF THEY APPEAR ONCE
# IN MOUSE GENOME

$|=1;				# flush print commands

my %hpeak=();
my %hchr=();

my $maxlen = shift @ARGV;
my $offset = shift @ARGV;

print STDERR "Reading tags... ";
while (<>){			# read one sequence at a time
    chomp;
    my ($id,$nmerge,$score,$chr,$start,$end,$lmerge,$mstart,$mend,@mergelist) = split/\t/; 
    if ($lmerge <= $maxlen) {
	$hpeak{$id}{score} = $score;
	$hpeak{$id}{chr} = $chr;
	$hpeak{$id}{start} = $start;
	$hpeak{$id}{end} = $end;
	$hchr{$chr}{$id}=1;
    }
}
my $npeak = keys %hpeak;
print STDERR "finished reading $npeak peaks with length smaller than $maxlen.\n";
print STDERR "Now expanding the region by +/- $offset and recording the sequences\n";

# # READ EACH CHROMOSOME SEQUENCE AND FIND THE TAGS
# # For reference: http://doc.bioperl.org/releases/bioperl-1.0/Bio/DB/Fasta.html

for my $chr ("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY","chrM") {

    $filename = $chr.".fa";
    print STDERR "Reading $chr... ";
    $db = Bio::DB::Fasta->new("../mouse_data/ucsc/chromFa/$filename");
    $obj = $db->get_Seq_by_id($chr);
    $chrlen = $obj->length;
    print STDERR "done reading $chrlen bases\n";

    for my $id (sort keys %{$hchr{$chr}}) {
	if (1) {
	    my $subseq = $obj->subseq(($hpeak{$id}{start}-$offset)=>($hpeak{$id}{end}+$offset));
	    print ">$id +/- $offset\n$subseq\n";
	}
    }
}

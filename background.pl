use Bio::Perl;
use Bio::DB::Fasta;

# READS A FILE WITH A SEQUENCE PER LINE and CHECK IF THEY APPEAR ONCE
# IN MOUSE GENOME

$|=1;				# flush print commands
my $pathname = "/home/alkan/c/gulayse/neurod2/mouse_data/ucsc/";
my %hprefix=(
    "ens"=>"ens",
    "refseq"=>"RefSeq_",
    );
my %hchr=();

my $database = shift @ARGV;
my $regtype = shift @ARGV;
my $regorder = shift @ARGV;
my $offset = shift @ARGV;
die "Run as: $0 <ens|refseq> <utr5|cds|utr3|intron> <0|1|2|..> <length-of-trailing-segments>\n"
    unless (defined $database && defined $regtype && defined $offset);

print STDERR "Generating background sequences by expanding ALL $regtype $regorder regions by +/- $offset bases\n";
print STDERR "Reading $regtype database...";
my $filename = $pathname . $hprefix{$database} . "Region.txt";

open(FP,$filename) or die "Can not open $filename\n";

my %hdatabase=();
my %hchr=();
my $cnt=0;
while(<FP>){
    chomp;
    ($chr,$start,$end,$txid,$sense,$reg,$order,$regionid) = split/\t/;
    if ( ($reg eq $regtype) && ($order eq $regorder) ) {
	$hdatabase{$regionid}{chr}=$chr;
	$hdatabase{$regionid}{start}=$start;
	$hdatabase{$regionid}{end}=$end;
	$hdatabase{$regionid}{sense}=$sense;
	$hdatabase{$regionid}{order}=$order;
	$hdatabase{$regionid}{length}=$end-$start;
	$hchr{$chr}{$regionid}=1;
	$cnt += 1;
    }
}
die "zero regions found !!!\n" unless ($cnt>0);

print STDERR "finished reading $cnt $reg_$regorder regions.\n";

# # READ EACH CHROMOSOME SEQUENCE AND FIND THE TAGS
# # For reference: http://doc.bioperl.org/releases/bioperl-1.0/Bio/DB/Fasta.html

for my $chr ("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY","chrM") {

    $filename = $chr.".fa";
    print STDERR "Reading $chr... ";
    $db = Bio::DB::Fasta->new("../mouse_data/ucsc/chromFa/$filename");
    $obj = $db->get_Seq_by_id($chr);
    $chrlen = $obj->length;
    print STDERR "done reading $chrlen bases. \n";
    for my $id (sort keys %{$hchr{$chr}}) {
	my $subseq = $obj->subseq(($hdatabase{$id}{start}-$offset)=>($hdatabase{$id}{end}+$offset));
	print ">$id +/- $offset\n$subseq\n";
    }
}

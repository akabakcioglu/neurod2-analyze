use Bio::Perl;
use Bio::DB::Fasta;

# READS A FILE WITH A SEQUENCE PER LINE and CHECK IF THEY APPEAR ONCE
# IN MOUSE GENOME

$|=1;				# flush print commands

%hmergeloc=();
%hchrlen=();
%hchroffset=();

$radius = shift @ARGV;

$nshift = 16;

# # READ EACH CHROMOSOME SEQUENCE AND FIND THE TAGS
# # For reference: http://doc.bioperl.org/releases/bioperl-1.0/Bio/DB/Fasta.html

# # READ EACH CHROMOSOME SEQUENCE AND FIND THE TAGS
# # For reference: http://doc.bioperl.org/releases/bioperl-1.0/Bio/DB/Fasta.html

$offset = 0;

for $chr ("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY","chrM") {

    $filename = $chr.".fa";
    print STDERR "Reading $chr... ";
    $db = Bio::DB::Fasta->new("../mouse_data/ucsc/chromFa/$filename");
    $obj = $db->get_Seq_by_id($chr);
    $hchrlen{$chr} = $obj->length;
    $hchroffset{$chr} = $offset; # all genome on a single axis
    $offset += $hchrlen{$chr};
    print STDERR "done reading $hchrlen{$chr} bases\n";
}


print STDERR "Reading tags... ";
while (<>){			# read one sequence at a time
    chomp;
    ($id,$nmerge,$score,$chr,$start,$end,$lmerge,$mstart,$mend,@mergelist) = split/\t/; # for merge_50
    $hmergeloc{$id} = $start + $hchroffset{$chr};
    $hmergeloc{$id} >>= $nshift;
}
$offset >>= $nshift;

$npeak = keys %hmergeloc;
print STDERR "finished reading $npeak peaks.\n";

@distr = ();
foreach $x (0..$offset) {
    push @distr,0;
}

foreach $id (sort {$hmergeloc{$a}<=>$hmergeloc{$b}} keys %hmergeloc) {
    foreach $x (($hmergeloc{$id}-$radius)..($hmergeloc{$id}+$radius)) {
	$distr[$x] += 1;
    }
}

foreach $x (0..$offset-1) {
    print STDOUT "$x\t$distr[$x]\n";
}
foreach $chr (sort keys %hchroffset) {
    $x = $hchroffset{$chr}>>$nshift;
    print STDERR  "$x\t0\n";
}



# Mart ja Mikk may 2017
# Database building script for PlasmidSeeker
# Takes in multifasta file with all plasmid sequences as input
# Database in the format of list files, each plasmid separate, with names.txt for keeping real names

use warnings;
use strict;
use Getopt::Long;
use FindBin;
use File::Basename;

# GLOBALS

my $fastafile;
my $glistmaker = "$FindBin::RealBin/GenomeTester4/glistmaker";
my $dir_location = "plasmid_db";
my $help;
my $version;
my $time = localtime;

# PARAMETERS

my $word = 20;
my $threads = 32;

# OPTIONS
GetOptions(
    'i=s' => \$fastafile,
    'd=s' => \$dir_location,
	'w=i' => \$word,
    'h' => \$help,
    'v' => \$version,
	't=s' => \$threads,
     ) or die printHelp()."\n";

if ($version){ die "PlasmidSeeker builder v1.1 (July 2020)\n"; }

# Check presence of options and files
die printHelp()."\n" if(!$dir_location || !$fastafile);

die "Can't find glistmaker at $glistmaker" if !-e $glistmaker;
##########
### SUBS
##########

#General help
sub printHelp {
	print "Usage: $0 -i <PLASMIDS FASTA FILE>\n";
	print "Options:
	-i\t Input fasta file with all plasmid sequences
	-d\t Database directory (default „plasmid_db“, will be created if does not exist)
	-t\t Number of threads used (default 32)
	-w\t K-mer length used (default 20)\n
	-h\t Print this help
	-v\t Print version of the program\n";
	return "";
}

##########
### MAIN
##########

###### Step 1 - extracts all plasmids as single fna files to a tmp directory
print STDERR "Getting all individual plasmid sequences and putting them to $dir_location\_fna\n";
open(FNA,'<',$fastafile) or die("Please specify plasmid fasta file!");
system "mkdir $dir_location\_fna" if !-d "$dir_location\_fna";

my $count = 1;
while(<FNA>) {
	if($_ =~ /^>/) {
			close FILE;
			open(FILE,'>',"$dir_location\_fna\/plasmid_$count\.fna");
			$count++;
	}
	print FILE $_;
}

###### Step 2 - converts each fna file to a list file (multithreaded)
print STDERR "Converting plasmid fna sequences to list files in $dir_location\n";
my @fna_files = glob("$dir_location\_fna/*");
system "mkdir $dir_location" if !-d $dir_location;

$count = 0;
my $cmd;
my $total = 0;
foreach (@fna_files) {
	$count++;
	$cmd .= "$glistmaker $_ -w $word -o $dir_location/".basename($_)." & ";
	#Multithread execution if thread limit reached
	if ($count == $threads) {
		$total = $total + $threads;
		print STDERR "Done: $total\r";
		system "$cmd wait";
		$count = 0;
		$cmd = "";
	}
}
# Check for unexecuted commands (last nodes if number of threads is not a multiple of node count)
if ($cmd ne "") {
	system "$cmd wait";
}


###### Step 3 - build names.txt file to give all plasmids real names - default taken from FASTA headers
print STDERR "Creating names.txt file in $dir_location\n";
open(FILE,'>',"$dir_location/names.txt");
$count = scalar(@fna_files);
print FILE "# Database: $dir_location\tPlasmids total: $count\tBuilt on: $time\tK-mer length: $word\n#\n"; # Print all info to header
foreach(@fna_files) {
	chomp;
	my $cmd = "head -n 1 $_";
	my $head = `$cmd`;
	my $length = `awk '/^>/{if (l!="") print l; l=0; next}{l+=length(\$0)}END{print l}' $_`;
	#$head = (split(/\|\s+/,$head))[1]; Split only in RefSeq format
	chomp $head;
	chomp $length;
	print FILE basename($_)."\t$head\t$length\n";
}
close FILE;

## FINAL - remove fna files and dir, leave only db lists
#system "rm $dir_location\_fna/*.fna";
#system "rmdir $dir_location\_fna";

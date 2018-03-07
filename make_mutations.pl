# For making random mutations to plasmid fna files
# USAGE: make_mutations.pl [plasmid fna] [mutations number]

use warnings;
use strict;

open (FILE,"<",$ARGV[0]);

my $read;
my $char;
my $flag = 0;
my $count;
my @nucl = ("A","C","G","T");

my @plasmid = <FILE>;
my $head = shift @plasmid;

my $fna;

foreach(@plasmid) {
	chomp;
	$fna .= $_;
}

my $total = length($fna);

close FILE;

#Mutation locations
my %mutations;

for(my $i=1;$i<=$ARGV[1];$i++) {
	$mutations{int(rand($total))} = 1;
}

foreach(keys %mutations) {
	print STDERR "mut: $_, ";
	my $nc = substr($fna,$_,1,);
	my $new = $nucl[rand @nucl];
	while($nc eq $new) { $new = $nucl[rand @nucl]; }
	print STDERR "$nc -> $new\n";
	substr($fna,$_,1,$new);
}

my $i = 0;
print $head;
while($i<$total) {
	my $line = substr($fna,$i,80);
	print "$line\n";
	$i+=80;
}
#print "$fna\n";


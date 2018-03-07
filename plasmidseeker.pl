# Mikk & Mart 2017 apr
# Plasmidseeker [WGS list file] [plasmid list dir] [plasmid fasta dir] [closest bacterium genome] [Results file name]
#

use warnings;
use strict;
use Getopt::Long;

# K-mer and other program paths
my $glistmaker = "GenomeTester4/glistmaker";
my $glistquery = "GenomeTester4/glistquery";
my $glistcompare = "GenomeTester4/glistcompare";
my $gdistribution = "GenomeTester4/gdistribution";
my $rtest = "testfunction.R";

# INPUT FILES AND NAMES
my $wgs_list; # Raw reads input
my $dir_location; # Plasmid db directory
my $bacteria; # Closest bacterium genome to isolate - for finding isolate coverage
my $outputfile;
my $help;
my $version;
my %plasmid_ids; # For keeping database plasmid names. Plasmid names file in db - names.txt
my %results; # Result hash - all high coverage plasmids
my $bacterial_cov; # Isolate median coverage

# GLOBAL PARAMETERS
my $cluster_prc = 80; # Percent used to cluster plasmids 
my $prc_limit = 80; # Percent of plasmid unique k-mers that have to be found from sample
my $threads = 32; # Num of threads
my $word; # K-mer length, get from the database
my $verbose = 0;
my $ponly = 0;
my $read_length; # Read length, get from data
my $bacterial_distr = "tmp_bacteria.txt"; # File of bacterial k-mer distribution
my $plasmid_distr = "tmp_plasmid.txt"; # File of plasmid k-mer distribution
my $multiple_tests = 1; # Bonferroni correction - number of tests done with high cov plasmids
my $distr_min = 2; # Minimum depth of k-mers to retain - lesser are cut in later tests

# OPTIONS
GetOptions(
    'i=s' => \$wgs_list,
    'd=s' => \$dir_location,
	'b=s' => \$bacteria,
    'h' => \$help,
    'v' => \$version,
    'o=s' => \$outputfile,
	't=s' => \$threads,
	'f=s' => \$prc_limit,
	'verbose' => \$verbose,
	'ponly' => \$ponly
    ) or die printHelp()."\n";

if ($version){ die "PlasmidSeeker v0.1 (21 Apr 2017)\n"; }

# Check presence of options and files
if(!$dir_location || !$wgs_list || (!$ponly && !$bacteria)) { die printHelp()."\n"; }
if (!$outputfile) { $outputfile = "plasmidseeker_result.txt"; }
if(!(qx/command -v Rscript/)) { die "Rscript not found - please make sure Rscript is available from PATH"; }

# Get read length
$read_length = getReadLen($wgs_list);

#############
# SUBS
#############

#General help
sub printHelp {
	print "Usage: $0 -d <PLASMID DB DIR> -i <SAMPLE.fastq> -b <CLOSEST BACTERIA TO ISOLATE> -o <OUTPUT FILE>\n";
	print "Options:
	-i\t Input fastq file
	-o\t Output file name (default plasmidseeker_result.txt)
	-d\t Path to plasmid database directory
	-b\t Closest bacteria to isolate genome fna
	-t\t Number of threads used (default 32)
	-f\t Minimum threshold F - at least this fraction of unique k-mers that has to be found for a plasmid (default 80)\n
	-h\t Print this help
	--verbose\t Print out more working process
	--ponly\t Assumes that reads contain only plasmid sequences (use for extracted plasmids)
	-v\t Print version of the program\n";
	return "";
}

#Gets input read length
sub getReadLen {
	my $file = shift;
	open (FNA,$file) or die "Cannot open raw reads file for getting read length\n";
	my $i = 1;
	my $typeChar;
	my $lenString;
	while(<FNA>){
		chomp;
		if (!$typeChar){
			$typeChar = substr($_, 0, 1);
			next;
		}

		if ($typeChar eq substr($_, 0, 1) || substr($_, 0, 1) eq "+"){
			last;
		}
		$lenString .= $_;

		$i++;
	}
	close FNA;
	
	return length($lenString);
}

#Finds median of a gdistribution list. INPUT - gdistribution output
sub find_median { 
	my @dist = @_;
	if(scalar @dist < 1) { return 0; }
	my $total = 0;
	my @counts;
	foreach(@dist) {
		my @arr = split(/\t/,$_);
		$total += $arr[1];
	}
	$total = $total/2;
	foreach(@dist) {
		my @arr = split(/\t/,$_);
		$total -= $arr[1];
		if($total<0) { return $arr[0]; }
	}
}

#Cuts off both ends of a list (default 1,2 and 3x median). INPUT - two list file names. First is WGS reads and the other is the plasmid/bacteria which coverage in the reads is to be found
# TODO file notation and put everything in one dir and remove afterwards
sub cut_list {
	my $wgs = shift;
	my $distr_name = shift; # Name of distr file generated
	my $list = shift;
	my $zero_kmers = defined $_[0] ? shift : 0; # Checking for kmers not found
	system "$glistcompare $wgs $list -i -r first -o tmp_list_cov"; # Find intersection
	system "$glistcompare tmp_list_cov_$word\_intrsec.list tmp_list_cov_$word\_intrsec.list -i -c $distr_min -o tmp_min_$distr_min"; # Remove 1-2 cov k-mers
	my @found = qx/$gdistribution tmp_list_cov_$word\_intrsec.list tmp_list_cov_$word\_intrsec.list/; # Find median
	print_distribution($distr_name,$zero_kmers,@found); # Print distribution for testing
	my $median = find_median(@found);
	print STDERR "Initial median $median\n";
	my $max_limit = int(3*$median);
	system "$glistcompare tmp_min_$distr_min\_$word\_intrsec.list tmp_min_$distr_min\_$word\_intrsec.list -i -c $max_limit -o tmp_longtail"; # Find part that is larger than max limit
	system "$glistcompare tmp_min_$distr_min\_$word\_intrsec.list tmp_longtail_$word\_intrsec.list -d -o tmp_final"; # Subtract longtail from main list for final cut list
	@found = qx/$gdistribution tmp_final_$word\_0_diff1.list tmp_final_$word\_0_diff1.list/; # Find final median
	$median = find_median(@found);
	print STDERR "FINAL median $median\n";
	return $median;
}

# Prints k-mer distribution to format suitable for R function
sub print_distribution {
	my $name = shift;
	my $zero_kmers = shift;
	chomp $name;
	print STDERR "Printing distribution of $name\n";
	open (my $fh,">",$name);
	print $fh "katvus;n\n";
	# Print k-mers not found in sample
	if($zero_kmers) { print $fh "0;$zero_kmers\n"; }
	foreach(@_) {
		my $line = $_; 
		$line =~ s/\t/\;/g;
		print $fh $line;
	}
	close $fh;
}

# Gets the testingfunction values for plasmid-bacteria (bacterial data must exist before can use this!). Returns testval, pval and koondus in hash
sub test_plasmid {
	print STDERR "Testing plasmid...\n";
	my @arr = qx/Rscript $rtest $bacterial_distr $plasmid_distr $read_length/;
	foreach(@arr) { chomp; }
	my ($teststat,$pval,$koond) = ((split(/\s+/,$arr[0]))[1],(split(/\s+/,$arr[1]))[1],(split(/\s+/,$arr[2]))[1]);
	my %values = ('teststat'=>$teststat,'pval'=>$pval,'koond'=>$koond);
	return %values;
}

#Parse glistcompare or glistquery output. Input - 
sub parse_output {
	my $type = shift;
	my @counts = @_;
	my %tmp;
	my $name;
	if($type eq 'compare') {
		foreach(@counts) {
			chomp;
			if($_ =~ /^1\t/) { 
				$name = (split(/\t/,$_))[1];
				next;
			}
			if($_ =~ /^NUnique/) { 
				my $count = (split(/\t/,$_))[1];
				$tmp{$name}{'S_unique'} = $count;
				next;
			}
			if($_ =~ /^NTotal/) { 
				my $count = (split(/\t/,$_))[1];
				$tmp{$name}{'S_total'} = $count;
				next;
			}
		}
	}
	if($type eq 'query') {
		foreach(@counts) {
			chomp;
			if($_ =~ /^Statistics of/) { 
				$name = (split(/\s+/,$_))[2];
				next;
			}
			if($_ =~ /^NUnique/) { 
				my $count = (split(/\t/,$_))[1];
				$tmp{$name}{'P_unique'} = $count;
				next;
			}
			if($_ =~ /^NTotal/) { 
				my $count = (split(/\t/,$_))[1];
				$tmp{$name}{'P_total'} = $count;
				next;
			}
		}
	}
	return %tmp;
}

####### END SUBS #########

################
# MAIN PROGRAM #
################

# Read in database and plasmid names
print STDERR "Loading database...\n";

opendir my $dir, $dir_location or die "Cannot open directory $!";
my @files = grep { $_ ne '.' && $_ ne '..' && $_ ne 'names.txt'} readdir $dir;
closedir $dir;
open(FILE,'<',"$dir_location/names.txt") or die "Plasmid names file $dir_location/names.txt missing!";
my @data = <FILE>;
my $header = shift(@data);
shift(@data);
chomp $header;
$word = (split(/length\:\s/,$header))[1];
if(!$word) { die "Cannot get k-mer length from database names.txt file!"; }
foreach(@data) {
	chomp;
	my @line = split(/\t/,$_);
	$plasmid_ids{$line[0]} = $line[1];
}
close FILE;

#Converting sample reads
print STDERR "Converting sample reads to k-mers...\n";
my $cmd = "$glistmaker $wgs_list -w $word -o tmp_sample";
system ($cmd);
$wgs_list = "tmp_sample_$word\.list";
system "$glistcompare tmp_sample_$word\.list tmp_sample_$word\.list -i -c $distr_min -o tmp_sample_min\_$distr_min"; # Remove 1 cov k-mers

#Bacteria genome coverage finding
if(!$ponly) {
	print STDERR "Finding coverage of bacterial isolate...\n";
	$cmd = "$glistmaker $bacteria -w $word -o tmp_closest";
	system ($cmd);
	$bacterial_cov = cut_list(($wgs_list,$bacterial_distr,"tmp_closest_$word\.list"));
	print STDERR "Bacteria median coverage is $bacterial_cov\n";
	# If coverage is too small for cutting list and other stuff...
	if($bacterial_cov<3) { die "Bacteria median coverage is too low (less than 3). Higher coverage dataset is needed or use flag --ponly"; }
}

##################################
# Multithreaded check all files
##################################

my %highcov;
my $count = 0;
my $total = 0;
my $cmd_cmp;
my $cmd_query;
my $total_plasmids = scalar(@files);
foreach(@files) {
	$count++;
	$cmd_cmp .= "$glistcompare tmp_sample_min\_$distr_min\_$word\_intrsec.list $dir_location/$_ -i -r first --count_only --print_operation &";
	$cmd_query .= "$glistquery $dir_location/$_ -stat &";
	
	#Multithread execution if thread limit reached
	if ($count == $threads) {
		$total += $threads;
		my @compare = qx/$cmd_cmp wait/;
		my @query = qx/$cmd_query wait/;
		my %out_compare = parse_output('compare',@compare);
		my %out_query = parse_output('query',@query);
		foreach(keys %out_compare) {
			my $prc = sprintf("%.2f",$out_compare{$_}{'S_unique'}/$out_query{$_}{'P_unique'}*100);
			
			#Over the F threshold - add to final results
			if($prc>$prc_limit) { $highcov{$_} = {
				Found=>$out_compare{$_}{'S_unique'},
				Sample_total=>$out_compare{$_}{'S_total'},
				Total=>$out_query{$_}{'P_total'},
				Total_unique=>$out_query{$_}{'P_unique'},
				Percent=>$prc,
			};
			}
			#print "$_ == found $out_compare{$_}{'S_unique'} of $out_query{$_}{'P_unique'}\t$prc\%\n";
		}
		print STDERR "Plasmids done: $total of $total_plasmids\r";
		$count = 0;
		$cmd_cmp = "";
		$cmd_query = "";
	}
}

# Check for unexecuted commands (last nodes if number of threads is not a multiple of node count)
if ($cmd_cmp ne "") {
	my @compare = qx/$cmd_cmp wait/;
	my @query = qx/$cmd_query wait/;
	my %out_compare = parse_output('compare',@compare);
	my %out_query = parse_output('query',@query);
	foreach(keys %out_compare) {
		my $prc = sprintf("%.2f",$out_compare{$_}{'S_unique'}/$out_query{$_}{'P_unique'}*100);
		
		#Over the F threshold - add to final results
		if($prc>$prc_limit) { $highcov{$_} = {
			Found=>$out_compare{$_}{'S_unique'},
			Sample_total=>$out_compare{$_}{'S_total'},
			Total=>$out_query{$_}{'P_total'},
			Total_unique=>$out_query{$_}{'P_unique'},
			Percent=>$prc,
		};
		}
	}
}

####### END MULTITHREADED PART ########

#check high cov plasmids
$multiple_tests = scalar keys %highcov;
foreach(keys %highcov) {
	# All other plasmid data
	my $re = qr/_$word\.list/;
	my $name = (split(/$re/,$_))[0];
	$name = (split(/\//,$name))[-1];
	my $plasmid_id = $plasmid_ids{$name};
	$plasmid_id = (split(/\|\s+/,$plasmid_id))[-1];
	$results{$_}{'Cov'} = sprintf("%.2f",$highcov{$_}{'Sample_total'}/$highcov{$_}{'Total'});
	chomp($highcov{$_}{'Found'});
	chomp($highcov{$_}{'Sample_total'});
	chomp($highcov{$_}{'Total'});
	chomp($highcov{$_}{'Total_unique'});
	$results{$_}{'Found_sample'} = $highcov{$_}{'Found'};
	$results{$_}{'Sample_total'} = $highcov{$_}{'Sample_total'};
	$results{$_}{'Plasmid_total'} = $highcov{$_}{'Total'};
	$results{$_}{'Plasmid_unique'} = $highcov{$_}{'Total_unique'};
	$results{$_}{'Percent'} = $highcov{$_}{'Percent'};
	$results{$_}{'Covered'} = int(($highcov{$_}{'Percent'}/100)*$highcov{$_}{'Total_unique'});
	$results{$_}{'PlasmidID'} = $plasmid_id;
	
	# Median-based coverage
	my $zero_kmers = $results{$_}{'Plasmid_unique'} - $results{$_}{'Found_sample'}; # Kmers not found
	$results{$_}{'Median'} = cut_list($wgs_list,$plasmid_distr,$_,$zero_kmers);
	
	if(!$ponly) {
		$results{$_}{'Cov_real'} = sprintf("%.2f",$results{$_}{'Median'}/$bacterial_cov);
	} else {
		$results{$_}{'Cov_real'} = $results{$_}{'Median'};
	}
	
	# Check percent found with cut list - to prevent false positives from erroneous k-mers (1-2x cov k-mers from other organisms)
	$cmd = "$glistquery tmp_final_$word\_0_diff1.list -stat";
	my %cutlist = parse_output('query',qx/$cmd/);
	my $cut_coverage = sprintf("%.2f",$cutlist{"tmp_final_$word\_0_diff1.list"}{'P_unique'}/$results{$_}{'Plasmid_unique'}*100);
	if($cut_coverage <= $prc_limit) { delete $results{$_}; next; }
	
	# Re-calculate found k-mers after subtracting low coverage k-mers
	$results{$_}{'Percent'} = $cut_coverage;
	$results{$_}{'Found_sample'} = $cutlist{"tmp_final_$word\_0_diff1.list"}{'P_unique'};
	$results{$_}{'Sample_total'} = $cutlist{"tmp_final_$word\_0_diff1.list"}{'P_total'};
	$results{$_}{'Covered'} = int(($results{$_}{'Percent'}/100)*$results{$_}{'Plasmid_unique'});
	
	# High-coverage plasmids included in the output - add statistical test
	if(!$ponly) {
		my %values = test_plasmid();
		$results{$_}{'Teststat'} = $values{'teststat'};
		$results{$_}{'Pval'} = $values{'pval'};
		$results{$_}{'Koond'} = $values{'koond'};
	} else {
		$results{$_}{'Teststat'} = 'na';
		$results{$_}{'Pval'} = 'na';
		$results{$_}{'Koond'} = 'na';
	}
}

#####################
# Clustering results
#####################

my %done; # Checked plasmids hash
my %cluster; # Cluster hash
$count = 1;

# Sorting by plasmid size, largest first
my @sorted_plasmids = sort {
    $results{$b}{'Plasmid_total'} <=> $results{$a}{'Plasmid_total'}
} keys %results;

my %high_pval;

# Filtering high p-value plasmids
my $corrected_pval = sprintf("%.6f",0.05/$multiple_tests);
foreach(keys %results) {
	if($results{$_}{'Pval'} eq 'na') { $high_pval{$_} = $results{$_}{'Pval'}; }
	else {
		my $pval = sprintf("%.6f",$results{$_}{'Pval'});
		if($pval > $corrected_pval) { $high_pval{$_} = $results{$_}{'Pval'}; }
	}
}

print STDERR "\nClustering results...\n";

foreach(@sorted_plasmids) {
	if(exists $done{$_} || exists $high_pval{$_}) { next; } # Skips plasmids that are already in cluster or high p-val
	#Start cluster
	$cluster{$count}{$_} = $results{$_}{'Covered'};
	$done{$_} = 1;
	
	foreach my $a (@sorted_plasmids) { #Search all results
		if(exists $done{$a} || exists $high_pval{$a}) { next; }
		
		my $cmd_isect = "$glistcompare $_ $a -i -r first --count_only";
		my @isect = qx/$cmd_isect/;
		my $shared = (split(/\t/,$isect[0]))[1];
		
		my $smaller = $results{$a}{'Plasmid_unique'};
		if($results{$a}{'Plasmid_unique'}>$results{$_}{'Plasmid_unique'}) { $smaller = $results{$_}{'Plasmid_unique'}; }
		if($shared/$smaller*100>$cluster_prc) {
			$cluster{$count}{$a} = $results{$a}{'Covered'};
			$done{$a} = 1;
		}
	}
	$count++;
}

###############
### OUTPUT
###############

if(!%cluster && !%results) { print STDERR "Nothing found...\n"; }

# Results header
open(my $fh,'>',$outputfile);
print $fh "# K-mers found\tTotal kmers\t%Kmers found(F)\tCopy number\tP-value\tPlasmid ID\tCoverage\tList file\n# MinP $prc_limit\n";
if(!$ponly) {
	print $fh "# Estimated bacteria isolate median coverage $bacterial_cov\n";
	print $fh "# Number of tests: $multiple_tests\tSignificant p-value (initial 0.05) with correction: $corrected_pval\n#\n";
}

#Clustered plasmids
foreach(my $i=1;$i<=$count;$i++) {
	if(!exists($cluster{$i})) { next; }
	print $fh "\tPLASMID CLUSTER $i\n";
	my %single_cluster = %{$cluster{$i}};
	foreach my $name (sort { $single_cluster{$b} <=> $single_cluster{$a} } keys %single_cluster) {
		print $fh "$results{$name}{'Covered'}\t$results{$name}{'Plasmid_unique'}\t$results{$name}{'Percent'}%\t$results{$name}{'Cov_real'}\t$results{$name}{'Pval'}\t$results{$name}{'PlasmidID'}\t$results{$name}{'Median'}\t$name\n";
	}
}

#High P-value plasmids
if(%high_pval) { print $fh "\n\tHIGH P-VALUE PLASMIDS\n"; }
foreach(keys %high_pval) {
	if(!exists $done{$_}) {
		print $fh "$results{$_}{'Covered'}\t$results{$_}{'Plasmid_unique'}\t$results{$_}{'Percent'}%\t$results{$_}{'Cov_real'}\t$results{$_}{'Pval'}\t$results{$_}{'PlasmidID'}\t$results{$_}{'Median'}\t$_\n";
	}
}

close $fh;
print STDERR "Done!\n";

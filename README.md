# PlasmidSeeker
A k-mer based program for the identification of known plasmids from bacterial whole genome sequencing reads

=========================

# USING PLASMIDSEEKER TOOL

## 1. INSTALLING AND BUILDING DATABASE

- Make sure you have [PERL](http://learn.perl.org/installing/unix_linux.html) and [R](https://cran.r-project.org/doc/manuals/r-release/R-admin.html#Obtaining-R) installed
- Put GenomeTester 4 binaries (gdistribution, glistcompare, glistquery, glistmaker) to a directory named "GenomeTester4", which should be directly under the main directory, which contains Testfunction.R, plasmidseeker.pl and database_builder.pl
- Put all plasmid FASTA files together into a single multi-fasta file (using UNIX cat command, for example)
- Use database_builder.pl to create the database or download our database with 8,514 plasmids (k=20) from our department server (http://bioinfo.ut.ee/plasmidseeker/) or from FigShare: https://figshare.com/s/5f7b924544839f7d6e59
- Approximate time with 8,514 Refseq plasmids, k=20 with 32 cores and 512GB RAM was 11 minutes.
- *For simplified installation and testing, follow the readme under "example" directory.*

command line example: *perl database_builder.pl -i [multi-FASTA file with all plasmids] -d [database directory name, will be created]*

### Database builder options:
- *-i* - Input fasta file with all plasmid sequences
- *-d* - Database directory, which will be created
- *-t* - Number of threads used (default 32)
- *-w* - K-mer length used (default 20)

## 2. DETECTING PLASMIDS

- Download closest bacterial reference genome to your isolate of interest (one possible source is NCBI Refseq: https://www.ncbi.nlm.nih.gov/refseq/)
- Use plasmidseeker.pl to detect plasmids from your isolate samples
- Approximate time using the 8,514 Refseq plasmids database, k=20 with 32 cores and 512GB RAM and 140 Mbp sample is less than 3 minutes.
- command line example: *perl plasmidseeker.pl -d [your database dir] -i [your WGS sample FASTQ file] -b [close reference bacterium FASTA file] -o [output file name]*

### PlasmidSeeker options:
- *-i* - Input FASTQ file location
- *-o* - Output file name (default is "plasmidseeker_result.txt")
- *-d* - Path to plasmid database directory
- *-b* - Closest reference bacterium FASTA file location
- *-t* - Number of threads used (default 32)
- *-f* - Minimum threshold F - at least this fraction of unique k-mers that has to be found for a plasmid (default 80)
- *--verbose* - Prints out more of the working process
- *--ponly* - Assumes that reads contain only plasmid sequences (use for extracted plasmids)

## 3. EXAMINING OUTPUT

- Reference plasmids which share more than 80% of k-mers are presented in a single cluster (example "CLUSTER 1") and ordered by the percentage of unique k-mers found. Likely, only one plasmid of each cluster is present in the sample.
- "HIGH P-VALUE PLASMIDS" refer to plasmids whose copy numbers are similar to the isolate. These plasmids may be integrated or false positives.
- Output example: "E_coli_wgs.txt" (E. coli WGS sample); "P_aeruginosa_plasmid.txt" (P. aeruginosa plasmid)

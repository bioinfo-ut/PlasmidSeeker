#!bin/bash
#Script for testing plasmid identification from E. coli WGS samples

chmod 755 GenomeTester4/glistquery
chmod 755 GenomeTester4/glistmaker
chmod 755 GenomeTester4/glistcompare
chmod 755 GenomeTester4/gdistribution
chmod 755 plasmidseeker.pl
chmod 755 database_builder.pl
chmod 755 testfunction.R

#Download PlasmidSeeker database and E. coli files

echo "Downloading PlasmidSeeker database..."
wget http://bioinfo.ut.ee/plasmidseeker/plasmidseeker_db_w20.tar.gz

echo "Unpacking database..."
tar -vzxf plasmidseeker_db_w20.tar.gz
rm plasmidseeker_db_w20.tar.gz

echo "Downloading E. coli data..."
wget http://bioinfo.ut.ee/plasmidseeker/e_coli_sakai_ref.fna
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR193/000/ERR1937840/ERR1937840.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR193/004/ERR1937914/ERR1937914.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR193/001/ERR1937841/ERR1937841.fastq.gz

echo "Unpacking files..."
gunzip ERR1937840.fastq.gz
gunzip ERR1937914.fastq.gz
gunzip ERR1937841.fastq.gz

#Use plamidseeker to identify E. coli samples

echo "Identifying E. coli sample EC1"
perl plasmidseeker.pl -d db_w20 -i ERR1937840.fastq -b e_coli_sakai_ref.fna -o EC_1_results.txt

echo "Identifying E. coli sample EC2"
perl plasmidseeker.pl -d db_w20 -i ERR1937914.fastq -b e_coli_sakai_ref.fna -o EC_2_results.txt

echo "Identifying E. coli sample EC3"
perl plasmidseeker.pl -d db_w20 -i ERR1937841.fastq -b e_coli_sakai_ref.fna -o EC_3_results.txt

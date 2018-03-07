# AN EXAMPLE OF DETECING PLASMIDS FROM <i>E. COLI</i> WGS SAMPLES   
In this example we are detecting known plasmids from <i>E. coli</i> WGS samples using a <i>k=20</i> database with 8514 known reference plasmids, as referenced in the PlasmidSeeker manuscript.
  
To do the analysis, you need to:
* a) download PlasmidSeeker repository containing bins, scripts and readme files from [Github](https://github.com/bioinfo-ut/PlasmidSeeker)  
* b) execute the bash script that downloads database and samples and executes plasmidseeker.    
    
Make sure you have enough space for storing these files. FASTA files that are used in this example are ca 1 GB and the database ca 8.8 GB.

Use following command lines to perform the example analysis ("plasmidseeker_ecoli_test.sh" downloads and unpacks database and FASTQ files and executes scripts):
```  
git clone https://github.com/bioinfo-ut/PlasmidSeeker/
bash plasmidseeker_ecoli_test.sh
```  
   
The 3 results files (EC_1_results.txt, EC_2_results.txt, EC_3_results.txt) contain the respective results of the <i>E. coli samples</i> mentioned in the manuscript (Table 2). 

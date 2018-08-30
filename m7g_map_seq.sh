# General m7G-MaP-seq data analysis pipeline (single end sequencing)

#### Requirements
# cutadapt (added to the path) https://cutadapt.readthedocs.io/en/stable/index.html
# samtools (added to the path) http://www.htslib.org/download/
# bowtie2 (added to the path) http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
# getFreq2000.R script (added to the path)
# R (https://www.r-project.org/)
# Fastq files from m7G-MaP-Seq experiment
# Fasta file with the sequences of the RNAs of interest

#### Optional (Use only if you have a barcode in your sequencing library)
# Preprocessing script from RNAprobR, which can be found here: https://github.com/lkie/RNAprobBash 
# (also contain guide for downloading and using RNAprobr including preprocessing script)
# Collapse script, which can be found here http://people.binf.ku.dk/jvinther/data/RNA-seq/collapse.sh

##############
# Preparation for mapping 
# Make directory for analysis
mkdir directory
cd /directory
mkdir data
cd /data

# Put fastq files in the data directory
# If necessary concatenate fastq files from the same index into one fastq file
for i in {1..6}
do
zcat "$i"_file1.fastq.gz "$i"_file2.fastq.gz > "$i".fastq.gz
done
wait

# To download example fastq files and E. coli rRNA fasta file use: 
wget -nH --cut-dirs=3 -r --no-parent --reject "index.html*" -e robots=off http://people.binf.ku.dk/jvinther/data/m7G-map-seq/data/

# Remove adapters (depends on the method used for library preparation, here standard Illumina adapter)
# For NextSeq sequencing use --nextseq-trim=20
for i in {1..6}
do
mkdir data/$i
cd /data/$i
zcat /data/"$i".fastq.gz | cutadapt -a AGATCGGAAGAGCACACGTCT --nextseq-trim=20 - 2> cutadapt.error | gzip > reads_trimmed.fastq.gz &
done
wait

##############
# For the example data a 7 base barcode is present in the adapter used for 3' cDNA ligation
# If you do not have a barcode in your sequencing library, you should skip the preprocessing step below and 
# do this command instead:
for i in {1..6}
do
cd /data/$i
mkdir output_dir
mv reads_trimmed.fastq.gz /data/$i/output_dir/read1.fastq.gz
done
wait

# OPTIONAL: Preprocessing for library with barcode
PATH=$PATH:/path/to/RNAprobr/scripts # set the path to the scripts
for i in {1..6} 
do
cd /data/"$i"
preprocessing.sh -b NNNNNNN -t 15 -1 reads_trimmed.fastq.gz 
gzip ./output_dir/read1.fastq
done
wait


##############
# Mapping (single end data)

#Make bowtie2 index file from Fasta file containing the sequences to be analysed
cd /data/fastafile
bowtie2-build Coli_rRNA.fa Coli_rRNA

# High sensitivity mapping with bowtie2 (single read)
# Local alignment with shortend seed that allows for a mismatch, increased number of seed extension attempts and reseeding.
for i in {1..6}
do
cd /data/"$i"
bowtie2 --local -N 1 -D 20 -R 3 -L 15 -x Coli_rRNA -U reads_trimmed.fastq.gz 2>bowtie2.error | gzip > mapped.sam.gz
done
wait

############
# OPTIONAL: if your library has barcode and was processed with the preprocessing script, the reads can be collapsed on
# barcodes to remove PCR duplicates from the analysis. 
# If you do not have a barcode in your sequencing library, you should skip the Collapse step below and 
# do this command instead:

for i in {1..6}
do
cd /data/$i
mv reads_mapped.sam.gz /data/$i/output_dir/mapped.sam.gz
done
wait

# OPTIONAL: remove barcodes containing N (the step below will cause these reads to be removed from analysis)
for i in {1..6}
do
cd /data/"$i"/output_dir
grep -P '^.*\t[^N]{7}' barcodes.txt > barcodes_filtered.txt
done
wait

# OPTIONAL: Collapse reads on the barcodes, script will remove reads that map to the same RNA position and have identical barcode
for i in {1..6}
do
cd /data/"$i"/output_dir
collapse.sh /data/"$i"/mapped.sam.gz barcodes_filtered.txt > mapped.sam.gz 
done
wait






##############
# Pileup of mapped reads using mpileup

# Sorting and making bams, necessary for input into mpileup
for i in {1..6}
do
cd /data/"$i"/output_dir
samtools view -u -S mapped.sam.gz|samtools sort > sorted.bam

# prepare a file listing the the path/filenames of the indexes relevant for the analysis
# Bam list example: index 1, 2, 3 are controls and index 4, 5, 6 are treated
mkdir data/output
cd data/output

echo "" > bam_file.txt
for i in {1..6}
do
echo /data/"$i"/output_dir/sorted.bam >> bam_file.txt
done
wait



#mpileup using bam list, use the same fasta file which were used for mapping
samtools mpileup -A -d 300000 -f /data/fastafile/Coli_rRNA.fa -b bam_file.txt > analysis.mpileup


##############
#Parsing mpileup file with getFreq Rscript
#Save the getFreq2000.R file to your system

###
# start R
R

source("path/getFreq2000.R")                           #read in the getFreq function
setwd("data/output")
getFreq("analysis.mpileup",                            #path to input mpileup file
                     CC = c(0,0,0,1,1,1),              #definition of the samples control=0 and treated=1 (relates to the order of samples in mpileup file)
                     CSVout = "path/analysis_out.txt", #output file
                     minFrac = 0.001,                  #minimum fraction of mutations for analysis to be performed for position
                     minCounts = 1,                    #minimum number of counts for analysis to be performed for position
                     pvalThres = -1,                   #sets p-val treshold for analysis to be performed for position (-1: all analysed)
                     nCores = 50,                      #number of computer cores used in analysis (depends on computer used)
                     nSites = 500,                     #number of positions analysed in parallel
                     maxSites = 500000)                #number of positions analysed (test the function with small number)


#The getFreq function will produce a tab delimited output file that can be read into R 
#For description of the columns in the output, see getFreq2000_output_description.xlsx file




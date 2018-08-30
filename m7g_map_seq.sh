# General m7G-MaP-seq data analysis pipeline (single end sequencing)

#### Requirements
# cutadapt (added to the path) https://cutadapt.readthedocs.io/en/stable/index.html
# samtools (added to the path) http://www.htslib.org/download/
# bowtie2 (added to the path) http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
# Fastq files from m7G-MaP-Seq experiment
# Fasta file with the sequences of the RNAs of interest

##############
# Preparation for mapping 
# Make directory for analysis
mkdir directory
cd /directory

# If necessary concatenate fastq files from the same index
cd PATH_TO_FASTQ
zcat index1_file1.fastq.gz index1_file2.fastq.gz > index1.fastq.gz

# For example fastq files and E. coli rRNA fasta file use: wget -r http://people.binf.ku.dk/jvinther/data/m7G-seq-map/data


# Remove adapters (depends on the method used for library preparation, here standard Illumina adapter)
# For NextSeq sequencing use --nextseq-trim=20
cutadapt -a AGATCGGAAGAGCACACGTCT --nextseq-trim=20 index1.fastq.gz 2> index1.cutadapt.error | gzip > index1_trimmed.fastq.gz

##############
# Mapping (single end data)

# Make bowtie2 index file from Fasta file containing the sequences to be analysed
bowtie2-build RNA_sequences.fa RNA_sequences

# High sensitivity mapping with bowtie2 (single read)
# Local alignment with shortend seed that allows for mismatches, reseeding
bowtie2 --local -N 1 -D 20 -R 3 -L 15 -x RNA_sequences -U index1_trimmed.fastq.gz 2>index1_bowtie2.error | gzip > index1_mapped.sam.gz

##############
# Pileup of mapped reads using mpileup

# Sorting and making bams, necessary for input into mpileup
samtools view -u -S index1_mapped.sam.gz|samtools sort > index1.bam

# prepare a file listing the the path/filenames of the indexes relevant for the analysis
# Bam list example: index 1, 2, 3 are controls and index 4, 5, 6 are treated
echo "" > bam_file.txt
for lab_number in 1 2 3 4 5 6
do
echo path/index"$lab_number".bam >> bam_file.txt
done
wait

#mpileup using bam list 
samtools mpileup -A -d 300000 -f RNA_sequences.fa -b bam_file.txt > analysis.mpileup

cd /binf-isilon/vintherlab/jvinther/070916/scratch/reanalysis150118/BacteriaPlant/bam_files/ara_wt_rrna


##############
#Parsing mpileup file with getFreq Rscript
#Save the getFreq2000.R file to your system
R

source("path/getFreq2000.R")                           #read in the getFreq function

getFreq("path/analysis.mpileup",                       #path to input mpileup file
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




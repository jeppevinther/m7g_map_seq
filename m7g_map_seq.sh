#General m7G-MaP-seq data analysis pipeline

##############
#Preparation for mapping (single end data)

#If necessary concatenate fastq files from the same index
zcat index1_file1.fastq.gz index1_file2.fastq.gz > index1.fastq.gz

#Remove adapters (depends on the method used for library preparation, here standard Illumina adapter)
#We have used cutadapt, but other tools such as Trimmomatic can also be used 
#For NextSeq sequencing use --nextseq-trim=20
cutadapt -a AGATCGGAAGAGCACACGTCT --nextseq-trim=20 index1.fastq.gz 2> index1.cutadapt.error | gzip > index1_trimmed.fastq.gz

##############
#Mapping (single end data)

#Make bowtie2 index file from Fasta file containing the sequences to be analysed
bowtie2-build RNA_sequences.fa RNA_sequences

#High sensitivity mapping with bowtie2 (single read) 
bowtie2 --local -N 1 -D 20 -R 3 -L 15 -x RNA_sequences -U index1_trimmed.fastq.gz 2>index1_bowtie2.error | gzip > index_mapped.sam.gz

##############
#Pileup of mapped reads







#

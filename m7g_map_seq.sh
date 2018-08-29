#General m7G-MaP-seq data analysis pipeline

##############
#Preparation for mapping

#If necessary concatenate fastq files from the same index
zcat file1.fastq.gz file2.fastq.gz

#
nice cutadapt -a AGATCGGAAGAGCACACGTCT --nextseq-trim=20 - 2>srna_m7g."$lab_number".R1.cutadapt.error | gzip > srna_m7g_"$lab_number"_R1.fastq.gz & 

nice zcat /seqdata/krogh/jvinther/170907/FASTQ/H5KY3BGX3_casava_1_mismatches/N127/*_S"$lab_number"_L00*_R1_001.fastq.gz* | nice cutadapt -a AGATCGGAAGAGCACACGTCT --nextseq-trim=20 - 2>srna_m7g."$lab_number".R1.cutadapt.error | gzip > srna_m7g_"$lab_number"_R1.fastq.gz & 


#Trim adapters and low quality sequences


#

#m7G-MaP-seq data analysis pipeline

##############
#Preparation for mapping

#If necessary concatenate fastq files from the same index


zcat file1.fastq.gz file2.fastq.gz | nice cutadapt -a AGATCGGAAGAGCACACGTCT --nextseq-trim=20 - 2>srna_m7g."$lab_number".R1.cutadapt.error | gzip > srna_m7g_"$lab_number"_R1.fastq.gz & 


#Trim adapters and low quality sequences


#

# Nextflow_bulk_RNAseq
Nextflow pipeline to trim, QC, align and quantify bulk RNA-seq fastq files

![image](https://github.com/user-attachments/assets/7a29336a-3cc9-4cb3-88b4-fd50bd70bd49)


# SARS-CoV-2 RNAseq scripts
This README file describes the Nextflow pipeline and R scripts used in order to process the human RNAseq datasets with respiratory virus infection in SARS-CoV-2.

## Table of Content:

1. Python script <filename: download_srr.py>: Python script used to download all SRR submissions (.fastq file/s) from NCBI for user entered SRP-ID

2. Nextflow pipeline <filename: fasterqdump_trimgalore_bowtie2_stringtie.nf>: 

      a. Process TRIM_GALORE: To trim low quality reads from downloaded *.fastq file/s

      b. Process BOWTIE2: To align trimmed mRNA reads from *_trimmed.fastq files to host genome using bowtie2 into *.sam file/s

      c. Process SAMTOOLS: To convert *.sam file/s to *.bam, *.aligned.bam, *.aligned.sorted.bam (bam file/s with reads aligned to host genome) and *.unaligned.bam, *.unaligned.sorted.bam (bam file/s with reads that dont aligned to host genome)

      d. Process STRINGTIE: To convert *.aligned.sorted.bam file/s to *.gtf file/s using host reference annotation file

      e. Processes FASTQC_pretrim and FASTQC_posttrim: To perform QC on reads before and after trimming

      f. Process MULTIQC: To combine FASTQC reports into a single report

      g. Process GENERATE_GTF_PATHS: To generate filenames and full-path/s to *.gtf file/s needed to convert *.gtf file/s into count matrix using prepDE.py

      h. Process prepDE: To convert *.gtf file/s generated by StringTie into count matrix
   
<img width="1077" alt="image" src="https://github.com/user-attachments/assets/741282f3-6864-4e46-9ebe-b50d4e1cdf5c">

4. R script <filename: DESeq2.R>: R script to perform Differential Gene Expression using DESeq2 and downstream analyses on the count matrix obtained from step 2h.

5. R script <filename: WGCNA.R>: R script to perform WGCNA on read count matrix obtained from step 3.





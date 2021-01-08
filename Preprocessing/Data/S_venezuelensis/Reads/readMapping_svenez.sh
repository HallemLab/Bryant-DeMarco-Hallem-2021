# This script checks the qualitiy of the fastq files and performs an alignment to the Strongyloides venezuelensis cDNA transcriptome reference with Kallisto.
# To run this 'shell script' you will need to open your terminal and navigate to the directory where this script resides on your computer.
# This should be the same directory where you fastq files and reference fasta file are found.
# Change permissions on your computer so that you can run a shell script by typing: 'chmod u+x readMapping.sh' (without the quotes) at the terminal prompt 
# Then type './readMapping_svenez.sh' (without the quotes) at the prompt.  
# This will begin the process of running each line of code in the shell script.

# first use fastqc to check the quality of the fastq files:
fastqc *.gz -t 14

# build index from the reference fasta file 
kallisto index -i strongyloides_venezuelensis.PRJEB530.WBPS14.mRNA_transcripts.index strongyloides_venezuelensis.PRJEB530.WBPS14.mRNA_transcripts.fa

# map reads to the indexed reference host transcriptome

# Parasitic Females: Biological Replicates 1-6, Technical replicate set 1
kallisto quant -i strongyloides_venezuelensis.PRJEB530.WBPS14.mRNA_transcripts.index -o DRR106346 -t 14 DRR106346_1.fastq.gz DRR106346_2.fastq.gz&> DRR106346.log
kallisto quant -i strongyloides_venezuelensis.PRJEB530.WBPS14.mRNA_transcripts.index -o DRR106348 -t 14 DRR106348_1.fastq.gz DRR106348_2.fastq.gz&> DRR106348.log
kallisto quant -i strongyloides_venezuelensis.PRJEB530.WBPS14.mRNA_transcripts.index -o DRR106350 -t 14 DRR106350_1.fastq.gz DRR106350_2.fastq.gz&> DRR106350.log
kallisto quant -i strongyloides_venezuelensis.PRJEB530.WBPS14.mRNA_transcripts.index -o DRR106352 -t 14 DRR106352_1.fastq.gz DRR106352_2.fastq.gz&> DRR106352.log
kallisto quant -i strongyloides_venezuelensis.PRJEB530.WBPS14.mRNA_transcripts.index -o DRR106354 -t 14 DRR106354_1.fastq.gz DRR106354_2.fastq.gz&> DRR106354.log
kallisto quant -i strongyloides_venezuelensis.PRJEB530.WBPS14.mRNA_transcripts.index -o DRR106356 -t 14 DRR106356_1.fastq.gz DRR106356_2.fastq.gz&> DRR106356.log

# Parasitic Females: Biological Replicates 1-6, Technical replicate set 2
kallisto quant -i strongyloides_venezuelensis.PRJEB530.WBPS14.mRNA_transcripts.index -o DRR106347 -t 14 DRR106347_1.fastq.gz DRR106347_2.fastq.gz&> DRR106347.log
kallisto quant -i strongyloides_venezuelensis.PRJEB530.WBPS14.mRNA_transcripts.index -o DRR106349 -t 14 DRR106349_1.fastq.gz DRR106349_2.fastq.gz&> DRR106349.log
kallisto quant -i strongyloides_venezuelensis.PRJEB530.WBPS14.mRNA_transcripts.index -o DRR106351 -t 14 DRR106351_1.fastq.gz DRR106351_2.fastq.gz&> DRR106351.log
kallisto quant -i strongyloides_venezuelensis.PRJEB530.WBPS14.mRNA_transcripts.index -o DRR106353 -t 14 DRR106353_1.fastq.gz DRR106353_2.fastq.gz&> DRR106353.log
kallisto quant -i strongyloides_venezuelensis.PRJEB530.WBPS14.mRNA_transcripts.index -o DRR106355 -t 14 DRR106355_1.fastq.gz DRR106355_2.fastq.gz&> DRR106355.log
kallisto quant -i strongyloides_venezuelensis.PRJEB530.WBPS14.mRNA_transcripts.index -o DRR106357 -t 14 DRR106357_1.fastq.gz DRR106357_2.fastq.gz&> DRR106357.log

# L1s: Technical Replicates 1-2
kallisto quant -i strongyloides_venezuelensis.PRJEB530.WBPS14.mRNA_transcripts.index -o DRR029433 -t 14 DRR029433_1.fastq.gz DRR029433_2.fastq.gz&> DRR029433.log
kallisto quant -i strongyloides_venezuelensis.PRJEB530.WBPS14.mRNA_transcripts.index -o DRR029434 -t 14 DRR029434_1.fastq.gz DRR029434_2.fastq.gz&> DRR029434.log

# iL3s: Technical Replicates 1-2
kallisto quant -i strongyloides_venezuelensis.PRJEB530.WBPS14.mRNA_transcripts.index -o DRR029435 -t 14 DRR029435_1.fastq.gz DRR029435_2.fastq.gz&> DRR029435.log
kallisto quant -i strongyloides_venezuelensis.PRJEB530.WBPS14.mRNA_transcripts.index -o DRR029436 -t 14 DRR029436_1.fastq.gz DRR029436_2.fastq.gz&> DRR029436.log

# L3s from lung: Technical Replicates 1-2
kallisto quant -i strongyloides_venezuelensis.PRJEB530.WBPS14.mRNA_transcripts.index -o DRR029437 -t 14 DRR029437_1.fastq.gz DRR029437_2.fastq.gz&> DRR029437.log
kallisto quant -i strongyloides_venezuelensis.PRJEB530.WBPS14.mRNA_transcripts.index -o DRR029438 -t 14 DRR029438_1.fastq.gz DRR029438_2.fastq.gz&> DRR029438.log

# Young Adult FLF: Technical Replicates 1-2
kallisto quant -i strongyloides_venezuelensis.PRJEB530.WBPS14.mRNA_transcripts.index -o DRR029439 -t 14 DRR029439_1.fastq.gz DRR029439_2.fastq.gz&> DRR029439.log
kallisto quant -i strongyloides_venezuelensis.PRJEB530.WBPS14.mRNA_transcripts.index -o DRR029440 -t 14 DRR029440_1.fastq.gz DRR029440_2.fastq.gz&> DRR029440.log

# Free-living Females: Technical Replicates 1-2 
kallisto quant -i strongyloides_venezuelensis.PRJEB530.WBPS14.mRNA_transcripts.index -o DRR029441 -t 14 DRR029441_1.fastq.gz DRR029441_2.fastq.gz&> DRR029441.log
kallisto quant -i strongyloides_venezuelensis.PRJEB530.WBPS14.mRNA_transcripts.index -o DRR029442 -t 14 DRR029442_1.fastq.gz DRR029442_2.fastq.gz&> DRR029442.log

# Activated iL3s, 1 day and 5 day (treat as Biological Replicates)
kallisto quant -i strongyloides_venezuelensis.PRJEB530.WBPS14.mRNA_transcripts.index -o DRR029443 -t 14 DRR029443_1.fastq.gz DRR029443_2.fastq.gz&> DRR029443.log
kallisto quant -i strongyloides_venezuelensis.PRJEB530.WBPS14.mRNA_transcripts.index -o DRR029444 -t 14 DRR029444_1.fastq.gz DRR029444_2.fastq.gz&> DRR029444.log

# Eggs: Technical Replicates 1-2 
kallisto quant -i strongyloides_venezuelensis.PRJEB530.WBPS14.mRNA_transcripts.index -o DRR029445 -t 14 DRR029445_1.fastq.gz DRR029445_2.fastq.gz&> DRR029445.log
kallisto quant -i strongyloides_venezuelensis.PRJEB530.WBPS14.mRNA_transcripts.index -o DRR029282 -t 14 DRR029282_1.fastq.gz DRR029282_2.fastq.gz&> DRR029282.log


# summarize fastqc and kallisto mapping results using MultiQC
multiqc -d . 
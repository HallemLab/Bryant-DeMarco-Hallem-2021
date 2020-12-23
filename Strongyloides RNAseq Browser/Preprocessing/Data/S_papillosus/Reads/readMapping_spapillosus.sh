# This script checks the qualitiy of the fastq files and performs an alignment to the Strongyloides papillosus cDNA transcriptome reference with Kallisto.
# To run this 'shell script' you will need to open your terminal and navigate to the directory where this script resides on your computer.
# This should be the same directory where you fastq files and reference fasta file are found.
# Change permissions on your computer so that you can run a shell script by typing: 'chmod u+x readMapping.sh' (without the quotes) at the terminal prompt 
# Then type './readMapping_svenez.sh' (without the quotes) at the prompt.  
# This will begin the process of running each line of code in the shell script.

# first use fastqc to check the quality of the fastq files:
fastqc *.gz -t 14

# build index from the reference fasta file 
kallisto index -i strongyloides_papillosus.PRJEB525.WBPS14.mRNA_transcripts.index strongyloides_papillosus.PRJEB525.WBPS14.mRNA_transcripts.fa

# map reads to the indexed reference host transcriptome

# L1/L2s: Biological Replicates 1-2
kallisto quant -i strongyloides_papillosus.PRJEB525.WBPS14.mRNA_transcripts.index -o ERR1466827 -t 14 ERR1466827_1.fastq.gz ERR1466827_2.fastq.gz&> ERR1466827.log
kallisto quant -i strongyloides_papillosus.PRJEB525.WBPS14.mRNA_transcripts.index -o ERR1466828 -t 14 ERR1466828_1.fastq.gz ERR1466828_2.fastq.gz&> ERR1466828.log

# Free-living Females: Biological Replicates 1-2
kallisto quant -i strongyloides_papillosus.PRJEB525.WBPS14.mRNA_transcripts.index -o ERR1466829 -t 14 ERR1466829_1.fastq.gz ERR1466829_2.fastq.gz&> ERR1466829.log
kallisto quant -i strongyloides_papillosus.PRJEB525.WBPS14.mRNA_transcripts.index -o ERR1466830 -t 14 ERR1466830_1.fastq.gz ERR1466830_2.fastq.gz&> ERR1466830.log

# Free-living Males: Biological Replicates 1-2
kallisto quant -i strongyloides_papillosus.PRJEB525.WBPS14.mRNA_transcripts.index -o ERR1466831 -t 14 ERR1466831_1.fastq.gz ERR1466831_2.fastq.gz&> ERR1466831.log
kallisto quant -i strongyloides_papillosus.PRJEB525.WBPS14.mRNA_transcripts.index -o ERR1466832 -t 14 ERR1466832_1.fastq.gz ERR1466832_2.fastq.gz&> ERR1466832.log

# iL3s: Biological Replicates 1-2
kallisto quant -i strongyloides_papillosus.PRJEB525.WBPS14.mRNA_transcripts.index -o ERR1466833 -t 14 ERR1466833_1.fastq.gz ERR1466833_2.fastq.gz&> ERR1466833.log
kallisto quant -i strongyloides_papillosus.PRJEB525.WBPS14.mRNA_transcripts.index -o ERR1466834 -t 14 ERR1466834_1.fastq.gz ERR1466834_2.fastq.gz&> ERR1466834.log

# Parasitic Females: biological Replicate 1
kallisto quant -i strongyloides_papillosus.PRJEB525.WBPS14.mRNA_transcripts.index -o ERR1466835 -t 14 ERR1466835_1.fastq.gz ERR1466835_2.fastq.gz&> ERR1466835.log

# Parasitic L1/L2s: Biological Replicate 1
kallisto quant -i strongyloides_papillosus.PRJEB525.WBPS14.mRNA_transcripts.index -o ERR1466836 -t 14 ERR1466836_1.fastq.gz ERR1466836_2.fastq.gz&> ERR1466836.log

# summarize fastqc and kallisto mapping results using MultiQC
multiqc -d . 
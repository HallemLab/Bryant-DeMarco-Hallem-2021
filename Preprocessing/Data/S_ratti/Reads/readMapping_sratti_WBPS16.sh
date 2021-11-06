# This script checks the qualitiy of the fastq files and performs an alignment to the Strongyloides ratti cDNA transcriptome reference with Kallisto.
# To run this 'shell script' you will need to open your terminal and navigate to the directory where this script resides on your computer.
# This should be the same directory where you fastq files and reference fasta file are found.
# Change permissions on your computer so that you can run a shell script by typing: 'chmod u+x readMapping_sratti_WBPS16.sh' (without the quotes) at the terminal prompt 
# Then type './readMapping_sratti_WBPS16.sh' (without the quotes) at the prompt.  
# This will begin the process of running each line of code in the shell script.

# first use fastqc to check the quality of the fastq files:
fastqc *.gz -t 14

# build index from the reference fasta file 
kallisto index -i strongyloides_ratti.PRJEB125.WBPS16.mRNA_transcripts.index strongyloides_ratti.PRJEB125.WBPS16.mRNA_transcripts.fa

# map reads to the indexed reference host transcriptome

# Parasitic Females: Biological Replicates A-C, Technical replicate set 1
kallisto quant -i strongyloides_ratti.PRJEB125.WBPS16.mRNA_transcripts.index -o ERR299168 -t 14 ERR299168_1.fastq.gz ERR299168_2.fastq.gz&> ERR299168.log
kallisto quant -i strongyloides_ratti.PRJEB125.WBPS16.mRNA_transcripts.index -o ERR299169 -t 14 ERR299169_1.fastq.gz ERR299169_2.fastq.gz&> ERR299169.log
kallisto quant -i strongyloides_ratti.PRJEB125.WBPS16.mRNA_transcripts.index -o ERR299170 -t 14 ERR299170_1.fastq.gz ERR299170_2.fastq.gz&> ERR299170.log

# Free-living Females: Biological Replicates A-C, Technical replicate set 1
kallisto quant -i strongyloides_ratti.PRJEB125.WBPS16.mRNA_transcripts.index -o ERR299171 -t 14 ERR299171_1.fastq.gz ERR299171_2.fastq.gz&> ERR299171.log
kallisto quant -i strongyloides_ratti.PRJEB125.WBPS16.mRNA_transcripts.index -o ERR299172 -t 14 ERR299172_1.fastq.gz ERR299172_2.fastq.gz&> ERR299172.log
kallisto quant -i strongyloides_ratti.PRJEB125.WBPS16.mRNA_transcripts.index -o ERR299173 -t 14 ERR299173_1.fastq.gz ERR299173_2.fastq.gz&> ERR299173.log

# Parasitic Females: Biological Replicates A-C, Technical replicate set 2
kallisto quant -i strongyloides_ratti.PRJEB125.WBPS16.mRNA_transcripts.index -o ERR299174 -t 14 ERR299174_1.fastq.gz ERR299174_2.fastq.gz&> ERR299174.log
kallisto quant -i strongyloides_ratti.PRJEB125.WBPS16.mRNA_transcripts.index -o ERR299175 -t 14 ERR299175_1.fastq.gz ERR299175_2.fastq.gz&> ERR299175.log
kallisto quant -i strongyloides_ratti.PRJEB125.WBPS16.mRNA_transcripts.index -o ERR299176 -t 14 ERR299176_1.fastq.gz ERR299176_2.fastq.gz&> ERR299176.log

# Free-living Females: Biological Replicates A-C, Technical replicate set 2
kallisto quant -i strongyloides_ratti.PRJEB125.WBPS16.mRNA_transcripts.index -o ERR299177 -t 14 ERR299177_1.fastq.gz ERR299177_2.fastq.gz&> ERR299177.log
kallisto quant -i strongyloides_ratti.PRJEB125.WBPS16.mRNA_transcripts.index -o ERR299178 -t 14 ERR299178_1.fastq.gz ERR299178_2.fastq.gz&> ERR299178.log
kallisto quant -i strongyloides_ratti.PRJEB125.WBPS16.mRNA_transcripts.index -o ERR299179 -t 14 ERR299179_1.fastq.gz ERR299179_2.fastq.gz&> ERR299179.log

# Free-living Males
kallisto quant -i strongyloides_ratti.PRJEB125.WBPS16.mRNA_transcripts.index -o ERR225783 -t 14 ERR225783_1.fastq.gz ERR225783_2.fastq.gz&> ERR225783.log

# iL3s
kallisto quant -i strongyloides_ratti.PRJEB125.WBPS16.mRNA_transcripts.index -o ERR225784 -t 14 ERR225784_1.fastq.gz ERR225784_2.fastq.gz&> ERR225784.log

# summarize fastqc and kallisto mapping results using MultiQC
multiqc -d . 

echo "Finished"
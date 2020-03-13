# Download and check the dataset
wget <source of data>
# Yoruba data as an example
samtools flagstat HGDP00944.alt_bwamem_GRCh38DH.20181023.Yoruba.cram

# Extract unmapped reads
samtools faidx ref.fa
samtools view -bt GRCh38_full_analysis_set_plus_decoy_hla.fa.fai -S -b -f 4 HGDP00944.alt_bwamem_GRCh38DH.20181023.Yoruba.cram > unmapped.bam

# Separate unmapped reads
# paired-end
samtools fastq -f 12 unmapped.bam -1 bothEndr1.fastq -2 bothEndr2.fastq
# one-end 
samtools fastq -f 68 -F 8 unmapped.bam > oneEndr1.fastq
samtools fastq -f 132 -F 8 unmapped.bam > oneEndr2.fastq

# assembly with MEGAHIT
megahit -1 bothEndr1.fastq -2 bothEndr2.fastq -m 0.5 -t 12 -o megahit_results_bothEnd
megahit -r oneEndr1.fastq -o megahit_results_r1
megahit -r oneEndr2.fastq -o megahit_results_r2

# remove contaminants with Kraken2
# minikraken2_v1-8GB is a pre-build library including Archaea, Bacteria and Virus from Kraken2 website
tar zxvf minikraken2_v1_8GB_201904_UPDATE.tgz
# remove contaminants from contigs
kraken2 --db minikraken2_v1_8GB megahit_results_bothEnd/final.contigs.fa > bothEnd_kraken2.txt
grep U bothEnd_kraken2.txt > unclassified.txt
awk ‘{print $2}’ unclassified.txt > unclassified_readID.txt
grep -w -A 1 -f unclassified_readID.txt megahit_results_bothEnd/final.contigs.fa --no-group-separator > final_contigs_bothEnd.fa

kraken2 --db minikraken2_v1_8GB megahit_results_r1/final.contigs.fa > oneEnd_r1_kraken2.t$
grep U oneEnd_r1_kraken2.txt > unclassified_r1.txt
awk ‘{print $2}’ unclassified_r1.txt > unclassified_r1_readID.txt
grep -w -A 1 -f unclassified_r1_readID.txt megahit_results_r1/final.contigs.fa --no-group-separator > final_contigs_r1.fa

kraken2 --db minikraken2_v1_8GB megahit_results_r2/final.contigs.fa > oneEnd_r2_kraken2.t$
grep U oneEnd_r2_kraken2.txt > unclassified_r2.txt
awk ‘{print $2}’ unclassified_r2.txt > unclassified_r2_readID.txt
grep -w -A 1 -f unclassified_r2_readID.txt megahit_results_r2/final.contigs.fa --no-group-separator > final_contigs_r2.fa


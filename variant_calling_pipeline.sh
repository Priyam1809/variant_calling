#!/bin/bash
# Variant Calling Pipeline Script
# Author: Priyam1809
# Date: 2025-11-18
#
# This script demonstrates a full variant calling workflow using common bioinformatics tools.
# It is designed for Linux systems. For Windows, use WSL or adapt commands accordingly.
#
# Placeholders are used for file paths and names. Replace them with your actual data.

set -e

# 1. Create directories for organization
mkdir -p samples/SRR____1 samples/SRR____2
mkdir -p fastqc_output ref trimmed_data fastqc_trim_output alignment vcf

echo "[1] Directories created for samples, outputs, and references."

# 2. Quality check of raw reads using FastQC
# Input: read1.fastq, read2.fastq (replace with your actual filenames)
# Output: HTML and summary files in fastqc_output/
fastqc raw_data/read1.fastq raw_data/read2.fastq -o fastqc_output

echo "[2] FastQC completed. Check fastqc_output/ for .html and summary files."

# 3. Download and index the reference genome (Homo sapiens)
# Download reference genome in .fasta format from NCBI or other source
# Place it as ref/homo_sapiens.fasta
# Indexing creates .bwt, .pac, .ann, .amb, .sa files for BWA
bwa index ref/homo_sapiens.fasta

echo "[3] Reference genome indexed with BWA."

# 4. Preprocess/trim raw reads to remove low-quality bases and adapters using fastp
# Output: trimmed FASTQ files in trimmed_data/
fastp -i raw_data/read1.fastq -I raw_data/read2.fastq -o trimmed_data/read1_trimmed.fastq -O trimmed_data/read2_trimmed.fastq --detect_adapter_for_pe

echo "[4] Read trimming completed. Trimmed files in trimmed_data/."

# 5. Quality check of trimmed reads
fastqc trimmed_data/read1_trimmed.fastq trimmed_data/read2_trimmed.fastq -o fastqc_trim_output

echo "[5] FastQC on trimmed reads completed. Check fastqc_trim_output/."

# 6. Align reads to the reference genome using BWA-MEM
# Output: SAM file in alignment/
bwa mem ref/homo_sapiens.fasta trimmed_data/read1_trimmed.fastq trimmed_data/read2_trimmed.fastq > alignment/read.sam

echo "[6] Alignment completed. SAM file in alignment/."

# 7. View the content of the SAM file (optional)
samtools view alignment/read.sam | head

echo "[7] Displayed first lines of SAM file."

# 8. Convert SAM to BAM
samtools view -b alignment/read.sam > alignment/read.bam

echo "[8] Converted SAM to BAM."

# 9. Sort BAM file
samtools sort alignment/read.bam -o alignment/read_sorted.bam

echo "[9] Sorted BAM file."

# 10. Index sorted BAM file
samtools index alignment/read_sorted.bam

echo "[10] BAM file indexed."

# 11. Generate mpileup file using bcftools
# Output: mpileup file in vcf/
bcftools mpileup -Ov -f ref/homo_sapiens.fasta alignment/read_sorted.bam > vcf/read_sorted.mpileup

echo "[11] Mpileup file generated."

# 12. Call variants with bcftools
# Output: VCF file in vcf/
bcftools call -c -v vcf/read_sorted.mpileup > vcf/read_sorted.vcf

echo "[12] Variant calling completed. VCF file in vcf/."

# 13. Count total number of variants and INDELs
variant_count=$(grep -v '^##' vcf/read_sorted.vcf | wc -l)
indel_count=$(grep -v '^##' vcf/read_sorted.vcf | grep 'INDEL' | wc -l)
echo "[13] Total variants: $variant_count"
echo "[13] Total INDELs: $indel_count"

# 14. Annotation with SnpEff
# Install SnpEff and Java if not present (requires sudo)
# sudo apt update
# sudo apt install -y openjdk-21-jre
# sudo apt install snpeff
# java -version
#
# Download SnpEff if not installed:
# wget https://snpeff.odsp.astrazeneca.com/versions/snpEff_latest_core.zip
# unzip snpEff_latest_core.zip
# cd snpEff
#
# Download the human GRCh38.86 database
# snpEff download GRCh38.86
#
# Annotate VCF
snpEff -v GRCh38.86 vcf/read_sorted.vcf > vcf/read_annotated.vcf

echo "[14] Annotation completed. Annotated VCF in vcf/."

# 15. Generate SnpEff HTML summary
snpEff -v -stats vcf/snpEff_summary.html GRCh38.86 vcf/read_sorted.vcf > vcf/read_annotated.vcf

echo "[15] SnpEff HTML summary generated at vcf/snpEff_summary.html."

# 16. (Optional) Custom reference annotation
# See script comments for custom reference genome annotation steps.

# End of pipeline

echo "Pipeline completed. Check output directories for results."

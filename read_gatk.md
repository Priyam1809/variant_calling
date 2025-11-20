# GATK Repository — Purpose and How to Use

## Why this repository exists

This repository provides a clear, reproducible set of commands and scripts for single-sample and cohort-based variant calling using GATK and common bioinformatics tools. It is intended to help researchers and bioinformaticians:
- Preprocess raw sequencing reads (QC and trimming).
- Align reads to a reference genome.
- Perform best-practice GATK post-processing (read groups, sorting, duplicate marking, BQSR).
- Call variants with GATK HaplotypeCaller.
- Filter or recalibrate variants (hard filters or VQSR).
- Produce final annotated VCFs ready for downstream analysis.

Why use this repo?
- Consolidates commands and explanations so you can run the pipeline step-by-step.
- Uses container-friendly instructions (Docker) for reproducibility.
- Includes placeholders and guidance for Windows users (via WSL or Docker).

---

## Key tools and short descriptions
- `fastqc` — FastQ quality control (per-base quality, adapters, duplication levels). Output: `.html` and `.zip` reports.
- `fastp` — Fast, all-in-one FASTQ preprocessor (adapter trimming, quality trimming, reports). Output: trimmed FASTQ, HTML/JSON report.
- `bwa` — Read aligner (BWA-MEM) to align reads to a reference. Output: SAM (or pipe to BAM).
- `samtools` — Manipulate SAM/BAM files (convert, sort, index).
- `gatk` — Broad Institute tools for post-processing and variant calling (AddOrReplaceReadGroups, SortSam, MarkDuplicates, BaseRecalibrator, ApplyBQSR, HaplotypeCaller, VariantFiltration, VariantRecalibrator, ApplyVQSR, SelectVariants).
- `bcftools` — VCF filtering/manipulation utilities.
- `snpEff` — Variant annotation.

---

## File layout (examples used in scripts)
- `raw_data/` — raw FASTQ files (place reads here)
- `trimmed_data/` — trimmed FASTQ outputs
- `fastqc_output/` — fastqc raw read reports
- `fastqc_trim_output/` — fastqc trimmed read reports
- `resource/` — reference fasta and known-sites VCFs (dbSNP, Mills, etc.)
- `alignment/` — SAM/BAM outputs and metrics
- `vcf/` — VCF outputs

---

## Step-by-step commands (Linux / WSL)
Below are the commands arranged in the standard order. Replace `PLACEHOLDER` with your file names and `REF` with your reference paths.

1) Preprocessing

```bash
# 1.1 Quality check of raw reads
fastqc raw_data/*.fastq* -o fastqc_output

# 1.2 Trim adapters and low-quality bases
fastp -i raw_data/READ1.fastq.gz -I raw_data/READ2.fastq.gz \
  -o trimmed_data/READ1_trim.fastq.gz -O trimmed_data/READ2_trim.fastq.gz \
  --detect_adapter_for_pe --thread 4

# 1.3 FastQC on trimmed reads
fastqc trimmed_data/*_trim.fastq* -o fastqc_trim_output
```

2) Reference indexing and alignment

```bash
# 2.1 Index reference for BWA (run once per reference)
bwa index resource/REF.fa

# 2.2 Align reads (BWA-MEM)
bwa mem resource/REF.fa trimmed_data/READ1_trim.fastq.gz trimmed_data/READ2_trim.fastq.gz > alignment/sample.sam

# 2.3 Convert to BAM, sort, and index (samtools recommended)
samtools view -b alignment/sample.sam -o alignment/sample.bam
samtools sort alignment/sample.bam -o alignment/sample_sorted.bam
samtools index alignment/sample_sorted.bam
```

3) GATK post-processing and variant calling (can run inside Docker container)

Docker (recommended for reproducibility): mount project dir and run GATK container:

```bash
# Pull GATK container
sudo docker pull broadinstitute/gatk:4.4.0.0

# Run container interactively (mount /absolute/path/to/project)
# inside container, your project is available at /gatk/data
sudo docker run -v /absolute/path/to/project:/gatk/data -it broadinstitute/gatk:4.4.0.0 /bin/bash
# Then inside container: cd /gatk/data and run gatk commands below
```

GATK commands (examples):

```bash
# Add or replace read groups (if input is BAM)
gatk AddOrReplaceReadGroups -I alignment/sample.bam -O alignment/sample_rg.bam \
  -RGID SAMPLE_ID -RGLB LIBRARY -RGPL ILLUMINA -RGPU UNIT -RGSM SAMPLE_NAME

# Sort (coordinate)
gatk SortSam -I alignment/sample_rg.bam -O alignment/sample_sorted.bam -SO coordinate

# Collect alignment metrics
gatk CollectAlignmentSummaryMetrics -R resource/REF.fa -I alignment/sample_sorted.bam -O alignment/alignment_metrics.txt

# Mark duplicates (optionally remove duplicates)
gatk MarkDuplicates -I alignment/sample_sorted.bam -O alignment/sample_sorted_dedup.bam -M alignment/marked_dup_metrics.txt --REMOVE_DUPLICATES true

# Index deduplicated BAM
samtools index alignment/sample_sorted_dedup.bam

# Base recalibration (BQSR) - needs known-sites VCFs (dbSNP, Mills)
gatk BaseRecalibrator -I alignment/sample_sorted_dedup.bam -R resource/REF.fa \
  --known-sites resource/Homo_sapiens_assembly38.known_indels.vcf.gz \
  --known-sites resource/Homo_sapiens_assembly38.dbsnp.vcf.gz -O alignment/recal_data.table

# Apply BQSR
gatk ApplyBQSR -I alignment/sample_sorted_dedup.bam -R resource/REF.fa \
  --bqsr-recal-file alignment/recal_data.table -O alignment/sample_sorted_dedup_recal.bam

# Variant calling
gatk HaplotypeCaller -R resource/REF.fa -I alignment/sample_sorted_dedup_recal.bam -O vcf/unfiltered.vcf \
  --dbsnp resource/Homo_sapiens_assembly38.dbsnp.vcf.gz --annotation AlleleFraction --annotation VariantType

# Hard-filtering (example thresholds) - adjust to your data
gatk VariantFiltration -V vcf/unfiltered.vcf \
  --filter-expression "QD < 2.0" --filter-name "QD2" \
  --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
  --filter-expression "SOR > 3.0" --filter-name "SOR3" \
  --filter-expression "FS > 60.0" --filter-name "FS60" \
  --filter-expression "MQ < 40.0" --filter-name "MQ40" -O vcf/hardfiltered.vcf

# Extract PASS variants using bcftools
bcftools filter -i 'FILTER=="PASS"' vcf/hardfiltered.vcf -o vcf/hardfilter_PASS.vcf
```

4) (Optional) VQSR — Variant Quality Score Recalibration
- VQSR requires high-quality resource files and generally works best on cohorts (many samples). If you have a single sample, hard-filtering is typically used instead.
- Example pattern (INDEL then SNP): use `gatk VariantRecalibrator` to build models and `gatk ApplyVQSR` to apply them.

5) Split variants into SNPs and INDELs

```bash
gatk SelectVariants -V vcf/final_PASS.vcf -R resource/REF.fa -O vcf/final_snps.vcf -select-type SNP
gatk SelectVariants -V vcf/final_PASS.vcf -R resource/REF.fa -O vcf/final_indels.vcf -select-type INDEL
```

6) Annotation with snpEff

```bash
# Download/prepare snpEff DB (one-time)
# java -jar snpEff.jar download GRCh38.86

# Annotate
snpEff -v GRCh38.86 vcf/final_PASS.vcf > vcf/final_PASS_annotated.vcf
```

---

## Windows usage

Recommended approach: Use WSL2 (Ubuntu) or Docker to run the exact Linux commands above. Examples:

- WSL2
  1. Install WSL2 and Ubuntu from Microsoft Store.
  2. From WSL bash, follow the Linux commands above.

- Docker (PowerShell example)

```powershell
# Run GATK container and mount C:\path\to\project
docker run -v C:\path\to\project:/gatk/data -it broadinstitute/gatk:4.4.0.0 /bin/bash
# Inside container:
cd /gatk/data
# Run the same gatk commands
```

Notes for native Windows: some bioinformatics tools have native Windows builds, but overall WSL/Docker offers the simplest route to run the Linux-first commands with minimal modification.

---

## Expected outputs at each major step
- FastQC: `.html` report per FASTQ (QC summary)
- fastp: trimmed FASTQ files and a JSON/HTML report
- BWA: SAM (or directly pipe to BAM)
- samtools/GATK: coordinate-sorted BAM, deduplicated BAM, BAM index
- GATK HaplotypeCaller: raw VCF (`unfiltered.vcf`)
- VariantFiltration / VQSR: filtered or recalibrated VCFs
- snpEff: annotated VCF with added annotation fields (effects, genes, consequence)

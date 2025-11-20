#!/bin/bash
# GATK Step-by-Step Interactive Pipeline
# Author: Priyam1809
# Date: 2025-11-20
#
# This script lists and runs the common GATK-based variant calling steps one by one.
# Each step includes a short explanation and expected outputs. Replace placeholders
# with your actual file paths before running.
#
# Usage (Linux/WSL):
#   1) Edit placeholders below (SEARCH FOR "PLACEHOLDER").
#   2) Make executable: chmod +x gatk_step_by_step.sh
#   3) Run interactively: ./gatk_step_by_step.sh
#
set -euo pipefail

### ======= PLACEHOLDERS: update these before running ======= ###
RAW_DIR="raw_data"                           # directory containing raw FASTQ files
TRIM_DIR="trimmed_data"                      # directory for trimmed FASTQ files
FASTQC_OUT="fastqc_output"                   # fastqc output directory
FASTQC_TRIM_OUT="fastqc_trim_output"        # fastqc output for trimmed data
REF_DIR="resource"                           # directory holding reference FASTA and known-sites
ALIGN_DIR="alignment"                        # alignment outputs
VCF_DIR="vcf"                                # VCF outputs
SAMPLE_PREFIX="gatk_demo"                    # prefix of your sample files
REF_FA="${REF_DIR}/Homo_sapiens_assembly38.fasta"   # Reference FASTA (hg38)
DBSNP_VCF="${REF_DIR}/Homo_sapiens_assembly38.dbsnp.vcf.gz"  # dbSNP VCF (gzipped)
KNOWN_INDELS="${REF_DIR}/Homo_sapiens_assembly38.known_indels.vcf.gz" # known indels

mkdir -p "$RAW_DIR" "$TRIM_DIR" "$FASTQC_OUT" "$FASTQC_TRIM_OUT" "$REF_DIR" "$ALIGN_DIR" "$VCF_DIR"

echo "Directories ensured: $RAW_DIR, $TRIM_DIR, $FASTQC_OUT, $FASTQC_TRIM_OUT, $REF_DIR, $ALIGN_DIR, $VCF_DIR"

pause(){
  read -rp "Press Enter to run this step (or Ctrl+C to abort) ..."
}

# STEP 1: Quality check of raw reads with FastQC
echo "\nSTEP 1: FastQC on raw reads"
cat <<'EXPL'
What it does:
 - FastQC inspects sequencing reads for per-base quality, GC content, adapter content, overrepresented sequences, etc.
Expected output:
 - One .html report and one .zip per FASTQ in $FASTQC_OUT
How to run:
 - Provide raw FASTQ files in $RAW_DIR; filenames can be .fastq or .fastq.gz
EXPL
pause
fastqc "$RAW_DIR"/*_R1*.fastq* "$RAW_DIR"/*_R2*.fastq* -o "$FASTQC_OUT" || echo "FastQC ended with non-zero status; inspect outputs"
ls -1 "$FASTQC_OUT"/*.html || true

# STEP 2: Trimming reads with fastp
echo "\nSTEP 2: fastp trimming"
cat <<'EXPL'
What it does:
 - fastp removes adapters, trims low-quality bases, and writes trimmed FASTQ files plus an HTML/JSON report.
Expected output:
 - Trimmed reads: ${TRIM_DIR}/${SAMPLE_PREFIX}_R1_trimmed.fastq.gz and _R2_trimmed.fastq.gz
 - fastp report: HTML and JSON in current directory (or TRIM_DIR if specified)
EXPL
pause
# Example: adjust filenames to match your raw data
fastp -i "$RAW_DIR/${SAMPLE_PREFIX}_R1.fastq.gz" -I "$RAW_DIR/${SAMPLE_PREFIX}_R2.fastq.gz" \
  -o "$TRIM_DIR/${SAMPLE_PREFIX}_R1_trimmed.fastq.gz" -O "$TRIM_DIR/${SAMPLE_PREFIX}_R2_trimmed.fastq.gz" \
  --detect_adapter_for_pe --thread 4 || echo "fastp may have returned non-zero status"
ls -1 "$TRIM_DIR"/*_trimmed.fastq* || true

# STEP 3: FastQC on trimmed reads
echo "\nSTEP 3: FastQC on trimmed reads"
cat <<'EXPL'
What it does:
 - Re-run FastQC to confirm trimming improved quality and removed adapters.
Expected output:
 - HTML and zip reports in $FASTQC_TRIM_OUT
EXPL
pause
fastqc "$TRIM_DIR"/*_trimmed.fastq* -o "$FASTQC_TRIM_OUT" || echo "FastQC (trimmed) returned non-zero status"
ls -1 "$FASTQC_TRIM_OUT"/*.html || true

# STEP 4: Index reference for BWA (one-time)
echo "\nSTEP 4: BWA index reference"
cat <<'EXPL'
What it does:
 - Creates index files required by BWA for alignment (.amb, .ann, .bwt, .pac, .sa)
Expected output:
 - index files in same folder as REF_FA
When to run:
 - Run once per reference FASTA (unless already indexed)
EXPL
pause
bwa index "$REF_FA" || echo "bwa index failed; check REF_FA path"
ls -1 "${REF_FA}"* || true

# STEP 5: Align trimmed reads with bwa mem
SAM_OUT="$ALIGN_DIR/${SAMPLE_PREFIX}.sam"
echo "\nSTEP 5: bwa mem alignment -> $SAM_OUT"
cat <<'EXPL'
What it does:
 - Aligns paired-end reads to the reference producing a SAM file.
Expected output:
 - SAM file at $SAM_OUT
EXPL
pause
bwa mem "$REF_FA" "$TRIM_DIR/${SAMPLE_PREFIX}_R1_trimmed.fastq.gz" "$TRIM_DIR/${SAMPLE_PREFIX}_R2_trimmed.fastq.gz" > "$SAM_OUT"
ls -lh "$SAM_OUT" || true

# STEP 6: Convert SAM to BAM and add read groups (GATK/Picard requires RG)
BAM_RG="$ALIGN_DIR/${SAMPLE_PREFIX}_rg.bam"
SORTED_BAM="$ALIGN_DIR/${SAMPLE_PREFIX}_sorted.bam"
DEDUP_BAM="$ALIGN_DIR/${SAMPLE_PREFIX}_sorted_dedup.bam"
RECAL_BAM="$ALIGN_DIR/${SAMPLE_PREFIX}_sorted_dedup_recal.bam"
UNFILTERED_VCF="$VCF_DIR/unfiltered.vcf"

echo "\nSTEP 6: Convert SAM->BAM and AddOrReplaceReadGroups"
cat <<'EXPL'
What it does:
 - Converts SAM to BAM (binary), then adds Read Group tags (RGID, RGLB, RGPL, RGPU, RGSM) required by GATK.
Expected output:
 - BAM with RG tags: $BAM_RG
EXPL
pause
samtools view -b "$SAM_OUT" -o "$ALIGN_DIR/${SAMPLE_PREFIX}.bam"
# Replace RG values below to match your sample (PLACEHOLDER SRR21388960)
gatk AddOrReplaceReadGroups -I "$ALIGN_DIR/${SAMPLE_PREFIX}.bam" -O "$BAM_RG" -RGID SRR21388960 -RGLB SRR21388960 -RGPL ILLUMINA -RGPU unit1 -RGSM SRR21388960 || echo "AddOrReplaceReadGroups issue"
ls -lh "$BAM_RG" || true

# STEP 7: Sort BAM by coordinate
echo "\nSTEP 7: Sort BAM (coordinate)"
cat <<'EXPL'
What it does:
 - Sorts alignments by genomic coordinate for downstream processing.
Expected output:
 - Sorted BAM at $SORTED_BAM
EXPL
pause
gatk SortSam -I "$BAM_RG" -O "$SORTED_BAM" -SO coordinate
ls -lh "$SORTED_BAM" || true

# STEP 8: Collect alignment metrics
echo "\nSTEP 8: Collect alignment metrics"
cat <<'EXPL'
What it does:
 - Collects metrics like total reads, mapped reads, insert sizes.
Expected output:
 - Text report: $ALIGN_DIR/alignment_metrics.txt
EXPL
pause
gatk CollectAlignmentSummaryMetrics -R "$REF_FA" -I "$SORTED_BAM" -O "$ALIGN_DIR/alignment_metrics.txt" || echo "CollectAlignmentSummaryMetrics may need Picard"
ls -lh "$ALIGN_DIR/alignment_metrics.txt" || true

# STEP 9: Mark duplicates
echo "\nSTEP 9: MarkDuplicates (and optionally remove)"
cat <<'EXPL'
What it does:
 - Marks (and can remove) PCR/optical duplicates. Produces metrics file.
Expected output:
 - Deduplicated BAM: $DEDUP_BAM
 - Metrics: $ALIGN_DIR/marked_dup_metrics.txt
EXPL
pause
gatk MarkDuplicates -I "$SORTED_BAM" -O "$DEDUP_BAM" -M "$ALIGN_DIR/marked_dup_metrics.txt" --REMOVE_DUPLICATES true
ls -lh "$DEDUP_BAM" "$ALIGN_DIR/marked_dup_metrics.txt" || true
samtools index "$DEDUP_BAM"

# STEP 10: Base Quality Score Recalibration (BQSR)
echo "\nSTEP 10: BaseRecalibrator -> ApplyBQSR"
cat <<'EXPL'
What it does:
 - BaseRecalibrator examines patterns of covariation and produces a recalibration table.
 - ApplyBQSR adjusts base quality scores using that table.
Expected output:
 - Recalibration table: $ALIGN_DIR/recal_data.table
 - Recalibrated BAM: $RECAL_BAM
Requirements:
 - Known-sites VCFs (dbSNP, Mills) present in $REF_DIR and indexed.
EXPL
pause
# Index known-sites VCFs if needed (indexing commands depend on file format: bgzip+tabix)
# Run BaseRecalibrator (may require large memory)
gatk BaseRecalibrator -I "$DEDUP_BAM" -R "$REF_FA" --known-sites "$KNOWN_INDELS" --known-sites "$DBSNP_VCF" -O "$ALIGN_DIR/recal_data.table"

# Apply recalibration
gatk ApplyBQSR -I "$DEDUP_BAM" -R "$REF_FA" --bqsr-recal-file "$ALIGN_DIR/recal_data.table" -O "$RECAL_BAM"
ls -lh "$RECAL_BAM" || true

# STEP 11: Variant calling with HaplotypeCaller
echo "\nSTEP 11: HaplotypeCaller (per-sample)"
cat <<'EXPL'
What it does:
 - Calls SNPs and INDELs to produce an unfiltered VCF.
Expected output:
 - VCF: $UNFILTERED_VCF
EXPL
pause
# Run HaplotypeCaller (may need --native-pair-hmm-threads or java options)
gatk HaplotypeCaller -R "$REF_FA" -I "$RECAL_BAM" -O "$UNFILTERED_VCF" --dbsnp "$DBSNP_VCF" --annotation AlleleFraction --annotation VariantType || echo "HaplotypeCaller completed"
ls -lh "$UNFILTERED_VCF" || true

# STEP 12: Hard-filter variants (example thresholds)
echo "\nSTEP 12: VariantFiltration (hard filters)"
HARD_FILTERED_VCF="$VCF_DIR/hardfiltered.vcf"
cat <<'EXPL'
What it does:
 - Applies site-level filters to flag low-quality calls.
 - For single-sample analysis hard-filtering is typical; VQSR is better for cohorts.
Expected output:
 - Filtered VCF: $HARD_FILTERED_VCF
EXPL
pause
gatk VariantFiltration -V "$UNFILTERED_VCF" \
  --filter-expression "QD < 2.0" --filter-name "QD2" \
  --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
  --filter-expression "SOR > 3.0" --filter-name "SOR3" \
  --filter-expression "FS > 60.0" --filter-name "FS60" \
  --filter-expression "MQ < 40.0" --filter-name "MQ40" -O "$HARD_FILTERED_VCF"
ls -lh "$HARD_FILTERED_VCF" || true

# Extract PASS variants with bcftools
echo "\nSTEP 13: Extract PASS variants"
PASS_VCF="$VCF_DIR/hardfilter_PASS.vcf"
cat <<'EXPL'
What it does:
 - Selects only variants that passed all filters (FILTER column == PASS).
Expected output:
 - PASS VCF: $PASS_VCF
EXPL
pause
bcftools filter -i 'FILTER="PASS"' "$HARD_FILTERED_VCF" -o "$PASS_VCF" || echo "bcftools filter may need installation"
ls -lh "$PASS_VCF" || true

# STEP 14: Variant Recalibration (VQSR) — recommended for large cohorts
echo "\nSTEP 14: VQSR (INDEL then SNP) — run inside GATK container or with sufficient resources"
cat <<'EXPL'
What it does:
 - Builds a statistical model to recalibrate variant quality scores using known truth sets.
 - Requires many samples for robust modeling; may fail for single-sample small datasets.
Expected outputs:
 - Recalibration files (.recal, .tranches) and a recalibrated VCF
EXPL
pause
cat <<'NOTE'
Example commands (run inside GATK docker or with proper gatk environment):

# INDEL recalibration
gatk --java-options "-Xmx4g -Xms4g" VariantRecalibrator \
  -V "$PASS_VCF" --trust-all-polymorphic \
  -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
  -mode INDEL \
  -resource:mills,known=false,training=true,truth=true,prior=12 $MILLS_INDEL \
  -resource:dbsnp,known=true,training=false,truth=false,prior=2 $DBSNP_VCF \
  -O "$VCF_DIR/hardfiltered_indel.recal" --tranches-file "$VCF_DIR/hardfiltered_indel.tranches"

gatk --java-options "-Xmx4g -Xms4g" ApplyVQSR -V "$PASS_VCF" --recal-file "$VCF_DIR/hardfiltered_indel.recal" --tranches-file "$VCF_DIR/hardfiltered_indel.tranches" --mode INDEL -O "$VCF_DIR/hardfiltered_indel_recal.vcf"

# SNP recalibration (similar pattern) and ApplyVQSR to produce final VCF
NOTE

echo "(VQSR commands printed above — run manually with correct resource file paths)"

# STEP 15: Split final PASS VCF into SNPs and INDELs (optional)
echo "\nSTEP 15: Split variants by type"
pause
cat <<'EXPL'
Commands:
 gatk SelectVariants -V <final_PASS.vcf> -R <ref.fa> -O snps_only.vcf -select-type SNP
 gatk SelectVariants -V <final_PASS.vcf> -R <ref.fa> -O indels_only.vcf -select-type INDEL
EXPL

# Final note
echo "\nAll scripted steps are complete. Inspect outputs in: $FASTQC_OUT, $FASTQC_TRIM_OUT, $ALIGN_DIR, $VCF_DIR"
echo "If you want, edit this script to match your filenames and run non-interactively by removing 'pause' calls."

exit 0

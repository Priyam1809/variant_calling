# Variant Calling Pipeline

## Purpose of This Repository

This repository provides a comprehensive, step-by-step workflow for variant calling using high-throughput sequencing data. It is designed for researchers and bioinformaticians who want to:
- Assess the quality of raw sequencing reads
- Preprocess and trim reads
- Align reads to a reference genome
- Call genetic variants (SNPs and INDELs)
- Annotate variants with biological information

The pipeline is suitable for human genome analysis but can be adapted for other organisms by changing the reference genome and annotation database.

## Included Commands and Tools
- `fastqc`: Quality control of raw and trimmed reads
- `fastp`: Read trimming and filtering
- `bwa`: Read alignment to reference genome
- `samtools`: SAM/BAM file manipulation
- `bcftools`: Variant calling
- `snpEff`: Variant annotation

## How to Use

### Linux (Recommended)
1. Clone the repository:
	```bash
	git clone https://github.com/Priyam1809/variant_calling.git
	cd variant_calling
	```
2. Make the pipeline script executable:
	```bash
	chmod +x variant_calling_pipeline.sh
	```
3. Edit the script to update file paths and sample names as needed.
4. Run the pipeline:
	```bash
	./variant_calling_pipeline.sh
	```
5. Review output files in the created directories (e.g., `fastqc_output/`, `alignment/`, `vcf/`).

### Windows
- **Recommended:** Use Windows Subsystem for Linux (WSL) to run the pipeline as above.
- Alternatively, install each tool for Windows and run the commands in PowerShell or Git Bash, adapting file paths as needed.
- Some tools (e.g., `bwa`, `samtools`, `bcftools`, `snpEff`) have Windows versions or can be run via WSL or Docker.

## Notes
- Replace all placeholder file paths (e.g., `raw_data/read1.fastq`) with your actual data files.
- Ensure all required tools are installed and available in your system's PATH.
- For large datasets, ensure you have sufficient disk space and memory.

## Output
- Quality reports (`.html`, `.zip`)
- Trimmed FASTQ files
- Alignment files (`.sam`, `.bam`)
- Variant files (`.vcf`)
- Annotated variant files

## License
This repository is open for academic and research use. Please cite the tools used according to their respective licenses.
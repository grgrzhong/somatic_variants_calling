# Somatic Variants Calling Pipeline

This pipeline performs somatic variant calling from paired-end DNA sequencing data (tumor-normal or tumor-only). It processes FASTQ files through quality control, alignment, somatic mutation calling, copy number variation analysis, and clinical annotation.

## Overview

The Somatic Variants Calling Pipeline includes the following main steps:
1. **Preprocessing & QC**: Quality control (FastQC) and adapter trimming with UMI extraction (fastp)
2. **BWA Alignment**: BWA-MEM alignment, UMI tagging, duplicate marking, and base quality recalibration
3. **Mutect2 Variant Calling**: GATK Mutect2 somatic variant calling with filtering and normalization
4. **CNV Analysis**: Allel specific copy number variation detection using FACETS
5. **Clinical Annotation**: Variant annotation using PCGR (Personal Cancer Genome Reporter)

## Installation

### Prerequisites

- Linux/Unix system with Singularity/Apptainer installed
- GNU Parallel for efficient sample processing
- Conda/Mamba for environment management

### Step 1: Install Conda/Mamba

```bash
# Download and install Miniforge (includes conda and mamba)
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash Miniforge3-Linux-x86_64.sh

# Follow the installation prompts and restart your terminal
source ~/.bashrc

# Verify installation
conda --version
```

### Step 2: Create Environment

```bash
# Create environment with required tools
conda create -n somatic_calling apptainer parallel
conda activate somatic_calling
```

### Step 3: Download and Setup the Pipeline

```bash
cd somatic_variants_calling

# Make the pipeline scripts executable
chmod +x somatic_calling_pipeline.sh
chmod +x run_pipeline.sh
```

## Pipeline Architecture

The pipeline consists of 5 main steps executed sequentially:

| Step | Script | Description | Key Tools |
|------|--------|-------------|-----------|
| 1 | `step_01_preprocessing.sh` | Quality control and adapter trimming with UMI extraction | fastp, FastQC |
| 2 | `step_02_bwa_alignment.sh` | BWA alignment, UMI tagging, sorting, duplicate marking, BQSR | BWA, GATK, SAMtools |
| 3 | `step_03_mutect2_call.sh` | Somatic variant calling, filtering, and normalization | GATK Mutect2, BCFtools |
| 4 | `step_04_cnv_facets.sh` | Copy number variation analysis | FACETS |
| 5 | `step_05_pcgr_annotation.sh` | Clinical variant annotation and reporting | PCGR, VEP |

## Configuration

### Main Configuration File: `conf/config.sh`

Edit the configuration file to set the following parameters:

- **`INPUT_DIR`**: Path to directory containing raw FASTQ files
- **`OUTPUT_DIR`**: Path to directory for output files  
- **`PARALLEL_JOBS`**: Number of samples to process in parallel (default: 1)
- **`REFERENCE_DIR`**: Path to reference files directory
- **`CONTAINER_DIR`**: Path to Singularity container files

### Required Reference Files

The pipeline requires several reference files in `REFERENCE_DIR`:

```
REFERENCE_DIR/
├── Gencode/
│   ├── gencode.hg38.v36.primary_assembly.fa       # Reference genome
│   └── annotation_protein_coding.bed              # Gene annotations
├── Exome/
│   └── xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged.bed  # Target intervals
├── Population_database/
│   ├── somatic-hg38_af-only-gnomad.hg38.vcf.gz   # Germline resource
│   └── dbSNP.vcf.gz                               # dbSNP variants
├── PON-Mutect/
│   ├── pon.vcf.gz                                 # Panel of normals
│   └── pon_pcgr.vcf.gz                           # PCGR-compatible PoN
├── PCGR_reference/                                # PCGR reference data
└── VEP_cache/                                     # VEP cache directory
```

## Usage

### Basic Usage

```bash
# Run the complete pipeline
./somatic_calling_pipeline.sh

# Or use the wrapper script
./run_pipeline.sh
```

### Custom Parameters

```bash
# Run with custom input/output directories
./somatic_calling_pipeline.sh \
  /path/to/input/fastq \
  /path/to/output \
  4  # parallel jobs
```

### Running Individual Steps

```bash
# Source configuration first
source conf/config.sh

# Run individual steps
bash scripts/workflow/step_01_preprocessing.sh
bash scripts/workflow/step_02_bwa_alignment.sh
bash scripts/workflow/step_03_mutect2_call.sh
bash scripts/workflow/step_04_cnv_facets.sh
bash scripts/workflow/step_05_pcgr.sh
```

## Input Requirements

### FASTQ File Naming Convention

The pipeline expects paired-end FASTQ files with the following naming pattern:
```
{sample_id}_1.fastq.gz  # Forward reads
{sample_id}_2.fastq.gz  # Reverse reads
```

### Sample Naming Convention

For tumor-normal pairs, use this naming scheme:
```
{PATIENT_ID}-T_1.fastq.gz    # Tumor sample
{PATIENT_ID}-T_2.fastq.gz
{PATIENT_ID}-N_1.fastq.gz    # Normal sample  
{PATIENT_ID}-N_2.fastq.gz
```

Example:
```
DFSP-001-T_1.fastq.gz
DFSP-001-T_2.fastq.gz
DFSP-001-N_1.fastq.gz
DFSP-001-N_2.fastq.gz
```

### Directory Structure Example

```
input_dir/
├── DFSP-001-T_1.fastq.gz
├── DFSP-001-T_2.fastq.gz
├── DFSP-001-N_1.fastq.gz
├── DFSP-001-N_2.fastq.gz
├── DFSP-002-T_1.fastq.gz
├── DFSP-002-T_2.fastq.gz
└── ...
```

## Output Structure

The pipeline generates the following output structure:

```
output_dir/
├── Fastq-trimmed/                    # Trimmed FASTQ files with UMI extraction
│   └── {sample}/
│       ├── {sample}_trimmed_1.fastq.gz
│       ├── {sample}_trimmed_2.fastq.gz
│       ├── {sample}.json             # fastp report
│       └── {sample}.html
├── FastQC-trimmed/                   # Quality control reports
│   └── {sample}/
│       └── *.html, *.zip
├── BAM/                              # Aligned and processed BAM files
│   └── {sample}/
│       ├── {sample}.recal.bam        # Final recalibrated BAM
│       ├── {sample}.recal.bam.bai    # BAM index
│       ├── {sample}.hsmetrics.txt    # Hybrid selection metrics
│       └── {sample}.stats.txt        # Alignment statistics
├── Mutect2/                          # Somatic variant calling results
│   └── {tumor_sample}/
│       ├── {tumor_sample}.final.vcf.gz           # Final filtered variants
│       ├── {tumor_sample}.mutect2.vcf.gz         # Raw Mutect2 calls
│       ├── {tumor_sample}.contamination.table    # Contamination estimates
│       └── *.log                                 # Process logs
├── CNV/                              # Copy number variation analysis
│   └── cnv_facets/
│       └── {tumor_sample}/
│           ├── {tumor_sample}.vcf.gz             # CNV calls in VCF format
│           ├── {tumor_sample}.tsv                # CNV calls in TSV format
│           └── {tumor_sample}.pcgr.tsv           # PCGR-compatible format
└── PCGR/                             # Clinical annotation results
    └── {tumor_sample}/
        ├── {tumor_sample}.pcgr.html              # Interactive HTML report
        ├── {tumor_sample}.pcgr.json              # Structured results
        ├── {tumor_sample}.pcgr.maf               # MAF format variants
        └── {tumor_sample}.pcgr.snvs_indels.tsv   # Annotated variants
```

## Key Features

### Tumor-Normal Mode
- Paired analysis for optimal somatic variant detection
- Contamination estimation and correction
- Improved filtering using matched normal samples

### Tumor-Only Mode  
- Single sample analysis when normal tissue unavailable
- Panel of Normals (PoN) filtering for artifact removal
- Population database filtering

### UMI Support
- Unique Molecular Identifier (UMI) extraction and tagging
- Improved duplicate marking accuracy
- Enhanced variant calling sensitivity

### Comprehensive Annotation
- Clinical variant interpretation using PCGR
- Tumor mutational burden (TMB) estimation
- Microsatellite instability (MSI) analysis
- Mutational signature analysis

## Quality Control

The pipeline includes multiple QC checkpoints:

1. **Raw data QC**: FastQC analysis of input FASTQ files
2. **Trimming QC**: Post-trimming quality assessment
3. **Alignment QC**: Hybrid selection metrics and alignment statistics
4. **Variant QC**: Multi-step filtering and normalization
5. **Contamination QC**: Cross-sample contamination detection

## Performance Considerations

- **Parallel Processing**: Use `PARALLEL_JOBS` parameter to process multiple samples simultaneously
- **Resource Requirements**: Each sample requires ~8-16 GB RAM and 8 CPU cores
- **Storage**: Plan for ~50-100 GB per sample for intermediate and final files
- **Runtime**: Complete analysis takes 4-8 hours per sample depending on coverage

## Troubleshooting

### Common Issues

1. **Missing reference files**: Ensure all required reference files are present and indexed
2. **Container access**: Verify Singularity containers are accessible and executable  
3. **Naming conventions**: Check that sample names follow the expected patterns
4. **Resource limits**: Monitor memory and CPU usage during execution

### Log Files

Check process-specific log files in each output subdirectory:
- `*.log` files contain detailed execution logs
- SLURM output files (if using cluster): `slurm/*.out` and `slurm/*.err`

## Citation

If you use this pipeline in your research, please cite the relevant tools:
- GATK Mutect2: Van der Auwera et al. (2013)
- BWA: Li & Durbin (2009) 
- FACETS: Shen & Seshan (2016)
- PCGR: Nakken et al. (2018)
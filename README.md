Somatic Variants Calling Pipeline

This pipeline performs somatic variant calling from paired-end DNA sequencing data (tumor-normal or tumor-only). It processes FASTQ files through quality control, alignment, somatic mutation calling, copy number variation analysis, and clinical annotation.

## Requirements

Install Conda and required software

```bash
## Download and install Miniforge (includes conda and mamba)
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash Miniforge3-Linux-x86_64.sh

## Follow the installation prompts and restart your terminal
source ~/.bashrc

## Verify installation
conda --version

## Create environment with required tools
conda install apptainer parallel

## Make the pipeline scripts executable
chmod +x somatic_calling_pipeline.sh
chmod +x run_pipeline_hpc.sh
```

## Input

Run the pipeline with Key parameters:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `INPUT_DIR` | Raw FASTQ files directory | `${PROJECT_DIR}/data/test_data/Raw` |
| `OUTPUT_DIR` | Results output directory | `${PROJECT_DIR}/results` |
| `PARALLEL_JOBS` | Parallel sample processing | `1` |
| `REFERENCE_DIR` | Reference files directory | `/lustre1/g/path_my/Reference` |
| `SOFTWARE_DIR` | Software directory | `/lustre1/g/path_my/Software` |
| `CONTAINER_DIR` | Singularity containers directory | `${PROJECT_DIR}/containers` |

The pipeline expects paired-end FASTQ files following this pattern:
```bash
{sample_id}_1.fastq.gz  # Forward reads (R1)
{sample_id}_2.fastq.gz  # Reverse reads (R2)
```


## Output

Example of outputs generated using this pipeline

```
├── Annovar
│   ├── DFSP-001-T
│   └── DFSP-030-T-P1
├── BAM
│   ├── DFSP-001-N
│   ├── DFSP-001-T
│   └── DFSP-030-T-P1
├── CNV_FACETS
│   └── DFSP-001-T
├── FastQC-trimmed
│   ├── DFSP-001-N
│   ├── DFSP-001-T
│   └── DFSP-030-T-P1
├── Fastq-trimmed
│   ├── DFSP-001-N
│   ├── DFSP-001-T
│   └── DFSP-030-T-P1
├── Mutect2
│   ├── DFSP-001-T
│   └── DFSP-030-T-P1
├── PCGR
│   ├── DFSP-001-T
│   └── DFSP-030-T-P1
├── Raw
│   ├── DFSP-001-N_1.fastq.gz
│   ├── DFSP-001-N_2.fastq.gz
│   ├── DFSP-001-T_1.fastq.gz
│   ├── DFSP-001-T_2.fastq.gz
│   ├── DFSP-030-T-P1_1.fastq.gz
│   └── DFSP-030-T-P1_2.fastq.gz
└── SV-call
    └── DFSP-001-T
```

## Quick start

```bash
## Custom parameters locally
./somatic_calling_pipeline.sh \
    /lustre1/g/path_my/pipeline/somatic_variants_calling/data/test_data/Raw \
    /lustre1/g/path_my/pipeline/somatic_variants_calling/data/test_data \
    2 \
    /lustre1/g/path_my/Reference \
    /lustre1/g/path_my/Software
```

```bash
## Submit to SLURM scheduler
./run_pipeline_hpc.sh
```
LineageXpress

LineageXpress is a probabilistic tool for predicting Mycobacterium tuberculosis complex (MTBC) lineages based on lineage-specific SNP markers.
It supports both sequential and parallel (multiprocessing) execution for efficient variant calling and lineage prediction on different HPC setups.

Features

SNP-based lineage classification using probabilistic scoring

Supports single-end, paired-end, BAM, or VCF inputs

Two execution modes:

Sequential version (sequential_version_lineage_predictor.py)

Parallel version (lineage_predictor_multiprocessing_pool_parallel.py)

LineageXpress/
├── README.MD
├── data
│   ├── h37rv.dict
│   ├── h37rv.fa
│   ├── h37rv.fa.*
│   ├── lineage_snp_updated_au13.tsv
│   └── targeted_modified_regions_au13.bed
├── environment_parallel.yml
├── environment_serial.yml
├── results/
├── sample_data
│   ├── SRR650226_1.fastq.gz
│   └── SRR650226_2.fastq.gz
├── sample_list.txt
└── scripts
    ├── parallel_version_lineageXpress.py
    └── sequential_version_lineageXpress.py

Create Conda Environments

# For serial version
conda env create -f environment/environment_serial.yml

# For parallel version
conda env create -f environment/environment_parallel.yml

Activate Environment
# Serial version
conda activate lineagexpress_serial

# Parallel version
conda activate lineagexpress_parallel

Usage

Sequential Version
Fastq input
python sequential_version_lineage_predictor.py 
  --fastq sample_list.txt \
  --output_dir results \

BAM input
python sequential_version_lineage_predictor.py 
  --fastq sample_list.txt \
  --output_dir results \

VCF input
python sequential_version_lineage_predictor.py 
  --fastq sample_list.txt \
  --output_dir results \

Parallel Version
python scripts/parallel_version_lineageXpress.py \
  --fastq sample_list.txt \
  --output_dir results \
  --n_jobs 1\
  --threads_per_tool 2

python scripts/parallel_version_lineageXpress.py \
  --bam sample_list.txt \
  --output_dir results \
  --n_jobs 1\
  --threads_per_tool 2

python scripts/parallel_version_lineageXpress.py \
  --vcf sample_list.txt \
  --output_dir results \
  --n_jobs 1 \
  --threads_per_tool 2


Quickstart Test (with provided sample data)

conda activate lineagexpress_serial 

python sequential_version_lineage_predictor.py 
  --fastq sample_list.txt \
  --output_dir results \
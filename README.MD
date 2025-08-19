<p align="center">
  <img src="assets/lineageXpress.png" alt="LineageXpress Logo" width="200"/>
</p>



LineageXpress is a probabilistic tool for predicting Mycobacterium tuberculosis complex (MTBC) lineages based on lineage-specific SNP markers.
It supports both sequential and parallel (multiprocessing) execution for efficient variant calling and lineage prediction on different setups, including HPC clusters.
Features

✨ Features

SNP-based lineage classification using probabilistic scoring

Supports single-end FASTQ, paired-end FASTQ, BAM, or VCF inputs

Two execution modes:

Sequential: sequential_version_lineage_predictor.py

Parallel: scripts/parallel_version_lineageXpress.py

Two execution modes:

Sequential version (sequential_version_lineage_predictor.py)

Parallel version (lineage_predictor_multiprocessing_pool_parallel.py)

Create Conda Environments

## For serial version
conda env create -f environment/environment_serial.yml

# For parallel version
conda env create -f environment/environment_parallel.yml

#Activate Environment
## Serial version
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





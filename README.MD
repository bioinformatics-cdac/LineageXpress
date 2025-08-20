<p align="center"
  <img src="assets/lineageXpress.png" alt="LineageXpress Logo" width="200"/
</p

# LineageXpress

**LineageXpress** is a probabilistic tool for predicting *Mycobacterium tuberculosis complex* (MTBC) lineages using lineage-specific SNP markers.  
It supports both **sequential** and **parallel (multiprocessing)** execution for efficient variant calling and lineage prediction on workstations and HPC clusters.

---

## ✨ Features
- SNP-based lineage classification using an interpretable probabilistic score
- Supports **single-end FASTQ**, **paired-end FASTQ**, **BAM**, or **VCF**
- Two execution modes:
  - **Sequential:** `sequential_version_lineageXpress.py`
  - **Parallel:** `parallel_version_lineageXpress.py`

---

## ⚙️ Installation (Conda)

```bash

git clone https://github.com/bioinformatics-cdac/LineageXpress.git
cd LineageXpress

# Serial (sequential) environment
conda env create -f environment_serial.yml
conda activate lineagexpress_serial

# Parallel (multiprocessing) environment
conda env create -f environment_parallel.yml
conda activate lineagexpress_parallel

```

## [!IMPORTANT] How to run this tool

 **Sequential**
 ```bash
 python scripts/sequential_version_lineageXpress.py \
   --fastq sample_list.txt \
   --output_dir results
 ```

 **Parallel**
 ```bash
 python scripts/parallel_version_lineageXpress.py \
   --fastq sample_list.txt \
   --output_dir results \
   --n_jobs 2 \
   --threads_per_tool 1
 ```
 Replace `--fastq` with `--bam` or `--vcf` to accept BAM/VCF inputs.

---

 🚀 Quickstart (with sample data)

 ```bash
 conda activate lineagexpress_serial
 python scripts/sequential_version_lineageXpress.py \
   --fastq sample_list.txt \
   --output_dir results
 ```
























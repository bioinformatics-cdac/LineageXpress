<p align="center">
  <img src="assets/lineageXpress.png" alt="LineageXpress Logo" width="200"/>
</p>

# LineageXpress

**LineageXpress** is a probabilistic tool for predicting *Mycobacterium tuberculosis complex* (MTBC) lineages using lineage-specific SNP markers.  
It supports both **sequential** and **parallel (multiprocessing)** execution for efficient variant calling and lineage prediction on workstations and HPC clusters.

---

## ✨ Features
- SNP-based lineage classification using an interpretable probabilistic score
- Supports **single-end FASTQ**, **paired-end FASTQ**, **BAM**, or **VCF**
- Two execution modes:
  - **Sequential:** `sequential_version_lineage_predictor.py`
  - **Parallel:** `scripts/parallel_version_lineageXpress.py`

---

## ⚙️ Installation (Conda)

```bash
# Serial (sequential) environment
conda env create -f environment/environment_serial.yml
conda activate lineagexpress_serial

# Parallel (multiprocessing) environment
conda env create -f environment/environment_parallel.yml
conda activate lineagexpress_parallel

<p align="center">
  <img src="https://raw.githubusercontent.com/bioinformatics-cdac/LineageXpress/master/assets/lineageXpress.png" 
       alt="LineageXpress Logo" width="200"/>
</p>

# LineageXpress

**LineageXpress** is a probabilistic tool for predicting *Mycobacterium tuberculosis complex* (MTBC) lineages using lineage-specific SNP markers.  
It supports both **sequential** and **parallel (multiprocessing)** execution for efficient variant calling and lineage prediction on workstations and HPC clusters. A simple **web GUI (Docker)** is also available.

---

## ✨ Features
- SNP-based lineage classification with an interpretable probability score
- Input types: **single-end FASTQ**, **paired-end FASTQ**, **BAM**, or **VCF**
- Execution modes:
  - **Sequential:** `sequential_version_lineageXpress.py`
  - **Parallel:** `parallel_version_lineageXpress.py`
- Optional **browser GUI** (Docker) with per-sample ZIP outputs

---

## ⚙️ Installation (Conda / CLI)

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

### Run (CLI)

**Sequential**
```bash
python scripts/sequential_version_lineageXpress.py --fastq sample_list.txt --output_dir results
```

**Parallel**
```bash
python scripts/parallel_version_lineageXpress.py --fastq sample_list.txt  --output_dir results   --n_jobs 2   --threads_per_tool 1
```

> Replace `--fastq` with `--bam` or `--vcf` to accept BAM/VCF inputs.

---

## 🚀 Quickstart (CLI, with sample data)

```bash
conda activate lineagexpress_serial
python scripts/sequential_version_lineageXpress.py --fastq sample_list.txt --output_dir results
```

---

## 🐳 Docker + Browser (GUI) — Beginner Friendly

The Docker image bundles the app and reference files. Just install Docker and follow these steps.

### 1) Pull the image (one-time)
```bash
docker pull bioinformaticscdac/lineagexpress:v1.0.0
```

### 2) Create a clean results folder on your machine
**Linux / WSL**
```bash
mkdir -p ~/lx_results
```
**Windows (PowerShell)**
```powershell
mkdir C:\lx_results
```

### 3) Start the web app
**Linux / WSL**
```bash
docker run --rm -it --name lxp-v1 -p 7861:7860 -v $HOME/lx_results:/app/results bioinformaticscdac/lineagexpress:v1.0.0
```

**Windows (PowerShell with Docker Desktop)**
```powershell
docker run --rm -it --name lxp-v1 -p 7861:7860  -v C:\lx_results:/app/results bioinformaticscdac/lineagexpress:v1.0.0
```

Open your browser at **http://localhost:7861**.

### 4) Use the GUI
1. Upload **one** of the following:
   - Paired FASTQ: select both **R1** and **R2**
   - **BAM** file
   - **VCF** file
2. Choose **Threads** (defaults are fine).
3. Click **Run Pipeline**.
4. When finished, click **Download All Results (ZIP)**.

```
> Tip: Don’t mount your entire home directory; use a dedicated folder like `~/lx_results` or `C:\lx_results`.

### 5) Stop the app
- In the web UI click **Close**, or press **Ctrl+C** in the terminal, or:
```bash
docker rm -f lxp-v1
```

### Run in background (optional)
```bash
docker run -d --name lxp-v1   --restart unless-stopped   -p 7861:7860   -v $HOME/lx_results:/app/results   bioinformaticscdac/lineagexpress:v1.0.0
```


```

---

📦 Outputs & File Layout

LineageXpress writes per-sample results to your chosen results/ directory (Docker: the host folder you mounted to /app/results).
Let SAMPLE_ID be the basename of your input (e.g., SRR650226).

1) Trimming & QC (Trim Galore + FastQC) (if trimming enabled)

Files:

SAMPLE_ID_1_val_1.fq.gz, SAMPLE_ID_2_val_2.fq.gz – trimmed reads

SAMPLE_ID_1_fastq.gz_trimming_report.txt, SAMPLE_ID_2_fastq.gz_trimming_report.txt – trimming summaries

(optional) SAMPLE_ID_1_fastqc.zip/html, SAMPLE_ID_2_fastqc.zip/html – FastQC reports

If you run the pipeline without trimming, these files won’t be produced.

2) Mapping to Reference (BWA + samtools)

Files:

SAMPLE_ID.sam – raw alignments (intermediate; large)

SAMPLE_ID.bam – BAM converted from SAM (intermediate)

SAMPLE_ID_mapped.bam – filtered/mapped reads (intermediate)

SAMPLE_ID_sorted.bam, SAMPLE_ID_sorted.bam.bai – final, indexed BAM

Logs:

logs/SAMPLE_ID.bwa.log – BWA/trim logs and stderr

3) Variant Calling (GATK HaplotypeCaller)

Files:

SAMPLE_ID.vcf – raw variants

SAMPLE_ID.vcf.idx – VCF index

(optional) SAMPLE_ID.filtered.vcf – variants after VariantFiltration (if enabled)

4) Lineage Prediction (SNP matching)

Files:

SAMPLE_ID_lineage_prediction_result.txt
(older builds: SAMPLE_ID_lineage_result.txt)

This text file includes:

Predicted lineage(s) and confidence

Mixed-lineage flag (if multiple high-confidence matches)

Per-lineage “matched/total SNPs” and probability

Total runtime

5) Per-sample ZIP (ready to share)

SAMPLE_ID_lineages_results.zip
Contains the key deliverables:

SAMPLE_ID_sorted.bam, SAMPLE_ID_sorted.bam.bai

SAMPLE_ID.vcf (and SAMPLE_ID.filtered.vcf if present)

SAMPLE_ID_lineage_prediction_result.txt

Intermediates like .sam, raw .bam, trimmed FASTQs, and FastQC zips can be deleted to save space after you’ve reviewed them.

## ✅ Troubleshooting

- **Browser doesn’t load:** confirm the port mapping shows in `docker ps` like `0.0.0.0:7861->7860/tcp`. Change host port if busy: `-p 9000:7860` → visit `http://localhost:9000`.
- **ZIP contains old files:** use a clean results folder (e.g., `~/lx_results`), not your whole home directory.
- **Slow performance:** increase the **Threads** slider if you have more CPU cores; using BAM/VCF input skips alignment.

---





import os
import subprocess
import logging
import argparse
import time
import re
from collections import defaultdict
from multiprocessing import Pool
from tqdm import tqdm

logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', level=logging.INFO)

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA_DIR = os.path.join(BASE_DIR, "data")
DEFAULT_DATA_DIR = "data"
DEFAULT_REF = os.path.join(DATA_DIR, "h37rv.fa")
DEFAULT_SNP_FILE = os.path.join(DATA_DIR, "lineage_snp_updated_au13.tsv")
DEFAULT_BED_FILE = os.path.join(DATA_DIR, "targeted_modified_regions_au13.bed")

def run_cmd(cmd):
    try:
        subprocess.run(cmd, check=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {cmd}\nSTDOUT: {e.stdout.decode()}\nSTDERR: {e.stderr.decode()}")
        return False
    except Exception as e:
        logging.error(f"An unexpected error occurred while running command: {cmd}\nError: {e}")
        return False
def load_lineage_snp_file(snp_file):
    lineage_snps = defaultdict(set)
    with open(snp_file, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                try:
                    parts = line.strip().split("\t")
                    if len(parts) == 3:
                        lineage, rvlocus, pos = parts
                    elif len(parts) == 2:
                        lineage, pos = parts
                        rvlocus = "NC_000962.3" # Default reference locus
                    else:
                        raise ValueError("Expected 2 or 3 columns")
                    lineage_snps[lineage].add((rvlocus, int(pos)))
                except ValueError as e:
                    logging.warning(f"Skipping malformed line: {line.strip()} â€” {e}")
    logging.info(f"Loaded SNP reference data for {len(lineage_snps)} lineages from {snp_file}")
    return lineage_snps

def parse_vcf_file(vcf_path):
    snps = set()
    try:
        with open(vcf_path, 'r') as vcf:
            for line in vcf:
                if not line.startswith('#'):
                    # Split only the first two parts to avoid issues with varying column numbers
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        chrom, pos = parts[0], parts[1]
                        try:
                            snps.add((chrom, int(pos)))
                        except ValueError:
                            logging.warning(f"Skipping VCF line with non-integer position: {line.strip()}")
                    else:
                        logging.warning(f"Skipping malformed VCF line (less than 2 columns): {line.strip()}")
    except FileNotFoundError:
        logging.error(f"VCF file not found: {vcf_path}")
        return set()
    return snps

def compare_snps(sample_id, sample_snps, lineage_references, output_file, elapsed_minutes):
    def calc_percentage(matched, total):
        return round((matched / total) * 100, 2) if total > 0 else 0.0

    results = {}
    matched_counts = {}
    for lineage, ref_snps in lineage_references.items():
        matched = len(ref_snps & sample_snps)
        total = len(ref_snps)
        probability = calc_percentage(matched, total)
        results[lineage] = probability
        matched_counts[lineage] = f"{matched}/{total}"

    formatted_results = {lin: f"{results.get(lin, 0.0):.2f}%" for lin in lineage_references}
    sorted_lineages = sorted(results.items(), key=lambda x: x[1], reverse=True)

    if not sorted_lineages or sorted_lineages[0][1] == 0.0:
        predicted_lineage = "Unpredictable: insufficient SNP evidence"
    else:
        top_pct = sorted_lineages[0][1]
        predicted_lineages = [lin for lin, pct in sorted_lineages if pct >= 90.0 or pct >= top_pct - 25.0]
        predicted_lineage = ";".join(sorted(set(predicted_lineages)))

    high_conf = [lin for lin, pct in results.items() if pct >= 90.0]
    mixed_details = ";".join(sorted(high_conf)) if len(high_conf) > 1 else "-"

    with open(output_file, 'w') as f:
        f.write(f"Sample              : {sample_id}\n")
        f.write(f"Predicted lineage  : {predicted_lineage}\n")
        f.write(f"Mixed lineage      : {mixed_details}\n\n")
        f.write(f"{'Lineage':<20}{'Probability':<15}{'Matched SNPs'}\n")
        f.write("-" * 50 + "\n")
        for lineage in sorted(lineage_references.keys()):
            prob = formatted_results.get(lineage, "0.00%")
            match = matched_counts.get(lineage, "0/0")
            f.write(f"{lineage:<20}{prob:<15}{match}\n")
        f.write("\n" + "=" * 60 + "\n")
        f.write(f"Total Runtime: {elapsed_minutes:.2f} minutes\n")

    logging.info(f"[{sample_id}] Lineage result saved to {output_file}")

def detect_fastq_files(fastq_dir, sample_base_id):
    patterns = [
        (f"{sample_base_id}_1.fastq", f"{sample_base_id}_2.fastq"),
        (f"{sample_base_id}_1.fastq.gz", f"{sample_base_id}_2.fastq.gz"),
        (f"{sample_base_id}_R1.fastq", f"{sample_base_id}_R2.fastq"),
        (f"{sample_base_id}_R1.fastq.gz", f"{sample_base_id}_R2.fastq.gz"),
        (f"{sample_base_id}.1.fastq", f"{sample_base_id}.2.fastq"),
        (f"{sample_base_id}.1.fastq.gz", f"{sample_base_id}.2.fastq.gz"),
    ]

    for r1_suffix, r2_suffix in patterns:
        r1_path = os.path.join(fastq_dir, r1_suffix)
        r2_path = os.path.join(fastq_dir, r2_suffix)
        if os.path.exists(r1_path) and os.path.exists(r2_path):
            return "paired", [r1_path, r2_path]

    for ext in [".fastq", ".fastq.gz"]:
        se_path = os.path.join(fastq_dir, sample_base_id + ext)
        if os.path.exists(se_path):
            return "single", [se_path]

    return "unknown", []
def process_fastq_single_sample(args):
    # Unpack tuple for multiprocessing compatibility
    (sample_entry, ref_genome, output_dir, lineage_references, bed_file, num_threads_per_tool) = args
    start_time = time.time()

    fastq_dir = os.path.dirname(sample_entry)
    base_filename = os.path.basename(sample_entry)
    
    sample_id = None
    match = re.match(r"(.*)(_R?[12]|\.[12])(\.fastq|\.fastq\.gz)$", base_filename)
    if match:
        sample_id = match.group(1)
    else:
        sample_id = os.path.splitext(base_filename)[0]
        if sample_id.endswith(".fastq") or sample_id.endswith(".fastq.gz"):
             sample_id = os.path.splitext(sample_id)[0]

    if not sample_id:
        logging.error(f"Could not determine sample ID from file: {sample_entry}")
        return

    sample_type, fastq_inputs = detect_fastq_files(fastq_dir, sample_id)

    if sample_type == "unknown":
        logging.error(f"[{sample_id}] Could not detect FASTQ files for sample: {sample_entry}")
        return

    logging.info(f"[{sample_id}] Starting processing of {sample_type} FASTQ files: {', '.join(fastq_inputs)}")

    logs_dir = os.path.join(output_dir, "logs")
    os.makedirs(logs_dir, exist_ok=True)

    sam_file = os.path.join(output_dir, f"{sample_id}.sam")
    bwa_log = os.path.join(logs_dir, f"{sample_id}.bwa.log")
    fastq_str = " ".join(fastq_inputs)
    bwa_cmd = (
        f"bwa mem -M -t {num_threads_per_tool} -R '@RG\\tID:{sample_id}\\tSM:{sample_id}\\tPL:ILLUMINA' "
        f"{ref_genome} {fastq_str} > {sam_file} 2> {bwa_log}"
    )
    if not run_cmd(bwa_cmd):
        logging.error(f"[{sample_id}] BWA command failed. Skipping further processing for this sample.")
        return

    bam_file = os.path.join(output_dir, f"{sample_id}.bam")
    mapped_bam = os.path.join(output_dir, f"{sample_id}_mapped.bam")
    sorted_bam = os.path.join(output_dir, f"{sample_id}_sorted.bam")

    if not run_cmd(f"samtools view -@ {num_threads_per_tool} -h -b -S -o {bam_file} {sam_file}"): return
    if not run_cmd(f"samtools view -@ {num_threads_per_tool} -h -b -F 260 -F 0x08 -F 0x400 -o {mapped_bam} {bam_file}"): return
    if not run_cmd(f"samtools sort -@ {num_threads_per_tool} {mapped_bam} -o {sorted_bam}"): return
    run_cmd(f"samtools index {sorted_bam}")

    vcf_file = os.path.join(output_dir, f"{sample_id}.vcf")
    gatk_cmd = f"gatk HaplotypeCaller -R {ref_genome} -I {sorted_bam} -O {vcf_file}"
    if bed_file:
        gatk_cmd += f" -L {bed_file}"
    if not run_cmd(gatk_cmd):
        logging.error(f"[{sample_id}] GATK HaplotypeCaller command failed. Skipping further processing for this sample.")
        return

    sample_snps = parse_vcf_file(vcf_file)
    output_file = os.path.join(output_dir, f"{sample_id}_lineage_result.txt")
    elapsed_minutes = (time.time() - start_time) / 60
    compare_snps(sample_id, sample_snps, lineage_references, output_file, elapsed_minutes)
    logging.info(f"[{sample_id}] Finished processing.")
def process_bam_single_sample(args):
    (bam_path, ref_genome, output_dir, lineage_references, bed_file, num_threads_per_tool) = args
    start_time = time.time()
    sample_id = os.path.splitext(os.path.basename(bam_path))[0]
    logging.info(f"[{sample_id}] Starting processing of BAM file.")

    sorted_bam = os.path.join(output_dir, f"{sample_id}.sorted.bam")
    if not run_cmd(f"samtools sort -@ {num_threads_per_tool} {bam_path} -o {sorted_bam}"): return
    if not run_cmd(f"samtools index {sorted_bam}"): return

    vcf_file = os.path.join(output_dir, f"{sample_id}.vcf")
    gatk_cmd = f"gatk HaplotypeCaller -R {ref_genome} -I {sorted_bam} -O {vcf_file}"
    if bed_file:
        gatk_cmd += f" -L {bed_file}"
    if not run_cmd(gatk_cmd):
        logging.error(f"[{sample_id}] GATK HaplotypeCaller command failed. Skipping further processing for this sample.")
        return

    sample_snps = parse_vcf_file(vcf_file)
    output_file = os.path.join(output_dir, f"{sample_id}_lineage_result.txt")
    elapsed_minutes = (time.time() - start_time) / 60
    compare_snps(sample_id, sample_snps, lineage_references, output_file, elapsed_minutes)
    logging.info(f"[{sample_id}] Finished processing.")

def process_vcf_single_sample(args):
    (vcf_path, output_dir, lineage_references) = args
    start_time = time.time()
    sample_id = os.path.splitext(os.path.basename(vcf_path))[0]
    logging.info(f"[{sample_id}] Starting processing of VCF file.")

    sample_snps = parse_vcf_file(vcf_path)
    output_file = os.path.join(output_dir, f"{sample_id}_lineage_result.txt")
    elapsed_minutes = (time.time() - start_time) / 60
    compare_snps(sample_id, sample_snps, lineage_references, output_file, elapsed_minutes)
    logging.info(f"[{sample_id}] Finished processing.")
    
def main():
    parser = argparse.ArgumentParser(
        description="LineageXpress: Predict MTBC lineage from FASTQ, BAM, or VCF"
    )

    # Mutually exclusive input types
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument("--fastq", help="Text file with list of FASTQ samples")
    group.add_argument("--bam", help="Text file with list of BAM samples")
    group.add_argument("--vcf", help="Text file with list of VCF samples")

    parser.add_argument("--output_dir", default="results",
                        help="Directory for outputs [default: results/]")
    parser.add_argument("--n_jobs", type=int, default=None,
                        help="Number of CPU cores to use [default: all available]")
    parser.add_argument("--threads_per_tool", type=int, default=1,
                        help="Threads per tool invocation [default: 1]")

    args = parser.parse_args()

    if not (args.fastq or args.bam or args.vcf):
        parser.print_help()
        print("\n Please provide one of: --fastq, --bam, or --vcf")
        return

    # Resolve reference files from data/
    ref_genome = DEFAULT_REF
    snp_file   = DEFAULT_SNP_FILE
    bed_file   = DEFAULT_BED_FILE

    # Load lineage SNP reference
    if not os.path.exists(snp_file):
        logging.error(f"SNP file missing: {snp_file}")
        return
    lineage_references = load_lineage_snp_file(snp_file)

    # Input file list
    input_file_list_path = args.fastq or args.bam or args.vcf
    if not os.path.exists(input_file_list_path):
        logging.error(f"Input file list not found: {input_file_list_path}")
        return

    with open(input_file_list_path, "r") as f:
        samples = [line.strip() for line in f if line.strip()]

    if not samples:
        logging.warning("No samples found in the input file list. Exiting.")
        return

    os.makedirs(args.output_dir, exist_ok=True)
    n_jobs = args.n_jobs or os.cpu_count()

    # Pick processing function
    if args.fastq:
        task_func = process_fastq_single_sample
        pool_args = [
            (sample, ref_genome, args.output_dir, lineage_references, bed_file, args.threads_per_tool)
            for sample in samples
        ]
        desc = "Processing FASTQ samples"
    elif args.bam:
        task_func = process_bam_single_sample
        pool_args = [
            (sample, ref_genome, args.output_dir, lineage_references, bed_file, args.threads_per_tool)
            for sample in samples
        ]
        desc = "Processing BAM samples"
    else:
        task_func = process_vcf_single_sample
        pool_args = [
            (sample, args.output_dir, lineage_references)
            for sample in samples
        ]
        desc = "Processing VCF samples"

    logging.info(f"Starting processing of {len(pool_args)} samples with {n_jobs} workers.")

    with Pool(processes=n_jobs) as pool:
        for _ in tqdm(pool.imap_unordered(task_func, pool_args), total=len(samples), desc=desc):
            pass

    logging.info("All samples processed.")


if __name__ == "__main__":
    main()
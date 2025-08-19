import os
import subprocess
import argparse
import logging
from collections import defaultdict
import time

logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', level=logging.INFO)

def run_cmd(cmd):
    try:
        subprocess.run(cmd, check=True, shell=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {cmd}\n{e}")
        return False
    return True

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
                        rvlocus = "NC_000962.3"
                    else:
                        raise ValueError("Expected 2 or 3 columns")
                    lineage_snps[lineage].add((rvlocus, int(pos)))
                except ValueError as e:
                    logging.warning(f"Skipping malformed line: {line.strip()} — {e}")
    logging.info(f"Loaded SNP reference data for {len(lineage_snps)} lineages from {snp_file}")
    return lineage_snps

def detect_sample_type(sample_path):
    dirname = os.path.dirname(sample_path)
    basename = os.path.basename(sample_path)
    sample_id = os.path.splitext(basename)[0]

    patterns = [
        (f"{sample_id}_1.fastq", f"{sample_id}_2.fastq"),
        (f"{sample_id}_1.fastq.gz", f"{sample_id}_2.fastq.gz"),
        (f"{sample_id}_R1.fastq", f"{sample_id}_R2.fastq"),
        (f"{sample_id}_R1.fastq.gz", f"{sample_id}_R2.fastq.gz"),
        (f"{sample_id}.1.fastq", f"{sample_id}.2.fastq"),
        (f"{sample_id}.1.fastq.gz", f"{sample_id}.2.fastq.gz"),
    ]

    for r1, r2 in patterns:
        r1_path = os.path.join(dirname, r1)
        r2_path = os.path.join(dirname, r2)
        if os.path.exists(r1_path) and os.path.exists(r2_path):
            return "paired", [r1_path, r2_path]

    for ext in [".fastq", ".fastq.gz"]:
        se_path = os.path.join(dirname, sample_id + ext)
        if os.path.exists(se_path):
            return "single", [se_path]

    return "unknown", []

def process_sample(sample_path, ref_genome, output_dir, lineage_references, bed_file, bam_override, vcf_override):
    def calc_percentage(matched, total):
        return round((matched / total) * 100, 2) if total > 0 else 0.0

    start_time = time.time()
    os.makedirs(output_dir, exist_ok=True)
    logs_dir = os.path.join(output_dir, "logs")
    os.makedirs(logs_dir, exist_ok=True)

    sample_id = os.path.splitext(os.path.basename(sample_path))[0]

    if vcf_override and os.path.exists(vcf_override):
        vcf_file = vcf_override
        logging.info(f"[{sample_id}] Using provided VCF file: {vcf_file}")
    else:
        if bam_override and os.path.exists(bam_override):
            sorted_bam = bam_override
            logging.info(f"[{sample_id}] Using provided BAM file: {sorted_bam}")
        else:
            sample_type, fastq_inputs = detect_sample_type(sample_path)
            if sample_type == "unknown":
                logging.error(f"[{sample_id}] Could not detect FASTQ files.")
                return

            sam_file = os.path.join(output_dir, f"{sample_id}.sam")
            bwa_log = os.path.join(logs_dir, f"{sample_id}.bwa.log")
            fastq_str = " ".join(fastq_inputs)
            bwa_cmd = (
                f"bwa mem -M -t 8 -R '@RG\\tID:{sample_id}\\tSM:{sample_id}\\tPL:ILLUMINA' "
                f"{ref_genome} {fastq_str} > {sam_file} 2> {bwa_log}"
            )
            logging.info(f"[{sample_id}] Running BWA MEM alignment")
            if not run_cmd(bwa_cmd) or not os.path.exists(sam_file):
                return

            bam_file = os.path.join(output_dir, f"{sample_id}.bam")
            sorted_bam = os.path.join(output_dir, f"{sample_id}.sorted.bam")

            if not run_cmd(f"samtools view -S -b {sam_file} -o {bam_file}"):
                return

            if not run_cmd(f"samtools sort {bam_file} -o {sorted_bam}"):
                return

            run_cmd(f"samtools index {sorted_bam}")

        vcf_file = os.path.join(output_dir, f"{sample_id}.vcf")
        if not run_cmd(
            f"gatk HaplotypeCaller -R {ref_genome} -I {sorted_bam} -O {vcf_file} -L {bed_file}"
        ):
            return

    sample_snps = set()
    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            if not line.startswith('#'):
                chrom, pos = line.strip().split('\t')[:2]
                sample_snps.add((chrom, int(pos)))

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

    elapsed_minutes = (time.time() - start_time) / 60
    output_file = os.path.join(output_dir, f"{sample_id}_lineage_result.txt")
    with open(output_file, 'w') as f:
        f.write(f"Sample             : {sample_id}\n")
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

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastq", required=True, help="Path to sample list file (one sample per line)")
    parser.add_argument("--ref_genome", default="data/h37rv.fa", help="Reference genome FASTA (default: data/h37rv.fa)")
    parser.add_argument("--output_dir", default="results", help="Output directory (default: results)")
    parser.add_argument("--snp_file", default="data/lineage_snps.txt", help="Lineage SNP reference file (default: data/lineage_snps.txt)")
    parser.add_argument("--bed_file", default="data/targets.bed", help="Target regions BED file (default: data/targets.bed)")
    parser.add_argument("--bam_file", default=None)
    parser.add_argument("--vcf_file", default=None)
    parser.add_argument("--threads", type=int, default=1, help="Threads for BWA (default: 1)")
    args = parser.parse_args()

    for path, label in [(args.ref_genome, "ref_genome"), (args.snp_file, "snp_file"), (args.bed_file, "bed_file")]:
        if not os.path.exists(path):
            raise FileNotFoundError(f"{label} not found at {path}. Place it in data/ or override the path.")
    lineage_references = load_lineage_snp_file(args.snp_file)
    with open(args.sample_file) as f:
        samples = [line.strip() for line in f if line.strip()]
        
   for sample_path in samples:
        process_sample(sample_path, args.ref_genome, args.output_dir,
                       lineage_references, args.bed_file,
                       bam_override=args.bam_file, vcf_override=args.vcf_file,threads=args.threads)

if __name__ == "__main__":
    main()


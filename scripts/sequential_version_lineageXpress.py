import os
import subprocess
import argparse
import logging
from collections import defaultdict
import time

logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', level=logging.INFO)

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA_DIR = os.path.join(BASE_DIR, "data")
DEFAULT_REF = os.path.join(DATA_DIR, "h37rv.fa")
DEFAULT_SNP_FILE = os.path.join(DATA_DIR, "lineage_snp_updated_au13.tsv")
DEFAULT_BED_FILE = os.path.join(DATA_DIR, "targeted_modified_regions_au13.bed")


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
def detect_sample_type(sample_entry):
    """
    Accepts:
      - file path (fastq/bam/vcf)
      - directory path (e.g., sample_data/SRR650226/)
      - bare ID (e.g., SRR650226) -> search under sample_data/<ID>/
    Returns: (kind, [paths]) where kind in {"paired","single","bam","vcf","unknown"}
    """
    # File given directly
    if os.path.isfile(sample_entry):
        name = sample_entry.lower()
        if name.endswith((".vcf", ".vcf.gz")):
            return "vcf", [sample_entry]
        if name.endswith(".bam"):
            return "bam", [sample_entry]
        if name.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz")):
            return "single", [sample_entry]

    # Directory or ID
    if os.path.isdir(sample_entry):
        base_dir = sample_entry
        sample_id = os.path.basename(os.path.normpath(sample_entry))
    else:
        base_dir = os.path.join("sample_data", sample_entry)
        sample_id = sample_entry

    # paired patterns
    patterns = [
        (f"{sample_id}_1.fastq",     f"{sample_id}_2.fastq"),
        (f"{sample_id}_1.fastq.gz",  f"{sample_id}_2.fastq.gz"),
        (f"{sample_id}_R1.fastq",    f"{sample_id}_R2.fastq"),
        (f"{sample_id}_R1.fastq.gz", f"{sample_id}_R2.fastq.gz"),
        (f"{sample_id}.1.fastq",     f"{sample_id}.2.fastq"),
        (f"{sample_id}.1.fastq.gz",  f"{sample_id}.2.fastq.gz"),
    ]
    for r1, r2 in patterns:
        r1_path = os.path.join(base_dir, r1)
        r2_path = os.path.join(base_dir, r2)
        if os.path.exists(r1_path) and os.path.exists(r2_path):
            return "paired", [r1_path, r2_path]

    # single-end
    for ext in (".fastq", ".fastq.gz", ".fq", ".fq.gz"):
        se = os.path.join(base_dir, sample_id + ext)
        if os.path.exists(se):
            return "single", [se]

    # BAM / VCF with sample_id basename
    bam = os.path.join(base_dir, f"{sample_id}.bam")
    if os.path.exists(bam):
        return "bam", [bam]
    for vext in (".vcf", ".vcf.gz"):
        v = os.path.join(base_dir, f"{sample_id}{vext}")
        if os.path.exists(v):
            return "vcf", [v]

    return "unknown", []
def process_sample(sample_path, ref_genome, output_dir, lineage_references,
                   bed_file, bam_override=None, vcf_override=None, threads=1):
    def pct(matched, total): return round((matched / total) * 100, 2) if total > 0 else 0.0

    start_time = time.time()
    os.makedirs(output_dir, exist_ok=True)
    logs_dir = os.path.join(output_dir, "logs")
    os.makedirs(logs_dir, exist_ok=True)

    sample_id = os.path.splitext(os.path.basename(str(sample_path)))[0]

    # Decide input source: VCF override > BAM override > detect
    vcf_file = None
    sorted_bam = None

    if vcf_override and os.path.exists(vcf_override):
        vcf_file = vcf_override
        logging.info(f"[{sample_id}] Using provided VCF: {vcf_file}")
    else:
        if bam_override and os.path.exists(bam_override):
            sorted_bam = bam_override
            logging.info(f"[{sample_id}] Using provided BAM: {sorted_bam}")
        else:
            kind, inputs = detect_sample_type(sample_path)
            if kind == "unknown":
                logging.error(f"[{sample_id}] Could not detect inputs under sample_data/{sample_id}")
                return
            if kind == "vcf":
                vcf_file = inputs[0]
                logging.info(f"[{sample_id}] Detected VCF: {vcf_file}")
            elif kind == "bam":
                sorted_bam = inputs[0]
                logging.info(f"[{sample_id}] Detected BAM: {sorted_bam}")
            else:
                # FASTQ -> align
                sam_file = os.path.join(output_dir, f"{sample_id}.sam")
                bwa_log = os.path.join(logs_dir, f"{sample_id}.bwa.log")
                fastq_str = " ".join(inputs)
                bwa_cmd = (
                    f"bwa mem -M -t {threads} -R '@RG\\tID:{sample_id}\\tSM:{sample_id}\\tPL:ILLUMINA' "
                    f"{ref_genome} {fastq_str} > {sam_file} 2> {bwa_log}"
                )
                logging.info(f"[{sample_id}] BWA MEM (threads={threads})")
                if not run_cmd(bwa_cmd) or not os.path.exists(sam_file):
                    return
                bam_file = os.path.join(output_dir, f"{sample_id}.bam")
                sorted_bam = os.path.join(output_dir, f"{sample_id}.sorted.bam")
                if not run_cmd(f"samtools view -S -b {sam_file} -o {bam_file}"): return
                if not run_cmd(f"samtools sort -@ {threads} {bam_file} -o {sorted_bam}"): return
                run_cmd(f"samtools index {sorted_bam}")

        # If we have a BAM but no VCF yet: call variants
        if sorted_bam and not vcf_file:
            vcf_file = os.path.join(output_dir, f"{sample_id}.vcf")
            if not run_cmd(f"gatk HaplotypeCaller -R {ref_genome} -I {sorted_bam} -O {vcf_file} -L {bed_file}"):
                return

    # Score
    sample_snps = set()
    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            if not line.startswith('#'):
                chrom, pos = line.strip().split('\t')[:2]
                try:
                    sample_snps.add((chrom, int(pos)))
                except ValueError:
                    pass

    results, matched_counts = {}, {}
    for lineage, ref_snps in lineage_references.items():
        matched = len(ref_snps & sample_snps)
        total = len(ref_snps)
        results[lineage] = pct(matched, total)
        matched_counts[lineage] = f"{matched}/{total}"

    formatted = {lin: f"{results.get(lin, 0.0):.2f}%" for lin in lineage_references}
    sorted_lineages = sorted(results.items(), key=lambda x: x[1], reverse=True)

    if not sorted_lineages or sorted_lineages[0][1] == 0.0:
        predicted = "Unpredictable: insufficient SNP evidence"
    else:
        top = sorted_lineages[0][1]
        predicted_set = [lin for lin, pc in sorted_lineages if pc >= 90.0 or pc >= top - 25.0]
        predicted = ";".join(sorted(set(predicted_set)))

    high_conf = [lin for lin, pc in results.items() if pc >= 90.0]
    mixed = ";".join(sorted(high_conf)) if len(high_conf) > 1 else "-"

    elapsed_minutes = (time.time() - start_time) / 60
    out = os.path.join(output_dir, f"{sample_id}_lineage_result.txt")
    with open(out, 'w') as f:
        f.write(f"Sample             : {sample_id}\n")
        f.write(f"Predicted lineage  : {predicted}\n")
        f.write(f"Mixed lineage      : {mixed}\n\n")
        f.write(f"{'Lineage':<20}{'Probability':<15}{'Matched SNPs'}\n")
        f.write("-" * 50 + "\n")
        for lineage in sorted(lineage_references.keys()):
            f.write(f"{lineage:<20}{formatted.get(lineage, '0.00%'):<15}{matched_counts.get(lineage, '0/0')}\n")
        f.write("\n" + "=" * 60 + "\n")
        f.write(f"Total Runtime: {elapsed_minutes:.2f} minutes\n")

    logging.info(f"[{sample_id}] Lineage result saved to {out}")
def main():
    parser = argparse.ArgumentParser(
        description="LineageXpress : Predict MTBC lineage from FASTQ, BAM, or VCF"
    )
    parser.add_argument("--sample_list", "--fastq", dest="sample_list", required=True,
                        help="Path to sample list file (one sample ID, dir, or file per line)")
    parser.add_argument("--ref_genome", default=DEFAULT_REF,
                        help=f"Reference genome FASTA (default: {DEFAULT_REF})")
    parser.add_argument("--output_dir", default="results",
                        help="Output directory (default: results)")
    parser.add_argument("--snp_file", default=DEFAULT_SNP_FILE,
                        help=f"Lineage SNP reference file (default: {DEFAULT_SNP_FILE})")
    parser.add_argument("--bed_file", default=DEFAULT_BED_FILE,
                        help=f"Target regions BED file (default: {DEFAULT_BED_FILE})")
    parser.add_argument("--bam_file", default=None,
                        help="Override: use this BAM for all samples")
    parser.add_argument("--vcf_file", default=None,
                        help="Override: use this VCF for all samples")
    parser.add_argument("--threads", type=int, default=1,
                        help="Threads for BWA / internal tools (default: 1)")
    args = parser.parse_args()

    # Preflight
    for path, label in [(args.ref_genome, "ref_genome"),
                        (args.snp_file, "snp_file"),
                        (args.bed_file, "bed_file")]:
        if path and not os.path.exists(path):
            raise FileNotFoundError(f"{label} not found at '{path}'. Place it in data/ or override the path.")

    lineage_refs = load_lineage_snp_file(args.snp_file)

    if not os.path.exists(args.sample_list):
        raise FileNotFoundError(f"sample_list file not found at '{args.sample_list}'")
    with open(args.sample_list, "r") as f:
        samples = [line.strip() for line in f if line.strip()]

    os.makedirs(args.output_dir, exist_ok=True)

    for sample in samples:
        process_sample(
            sample_path=sample,
            ref_genome=args.ref_genome,
            output_dir=args.output_dir,
            lineage_references=lineage_refs,
            bed_file=args.bed_file,
            bam_override=args.bam_file,
            vcf_override=args.vcf_file,
            threads=args.threads,
        )

if __name__ == "__main__":
    main()




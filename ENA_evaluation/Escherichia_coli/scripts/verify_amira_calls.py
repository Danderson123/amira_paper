import subprocess
import os
import argparse
import pysam
import statistics
import pandas as pd
import json
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument("--reads", dest="reads", required=True)
parser.add_argument("--amira-output", dest="amira_output", required=True)
parser.add_argument("--output", dest="output", required=True)
parser.add_argument("--cores", dest="cores", required=True)
args = parser.parse_args()

def map_reads_to_fasta(genes_file, nanopore, output_file, cores):
    # map the reads to the genes
    cmd = f"minimap2 -a -x map-ont -t {cores} --eqx --MD {genes_file} {nanopore} > {output_file}"
    output_bam = output_file.replace(".sam", ".bam")
    if not os.path.exists(output_file):
        subprocess.run(cmd, shell=True, check=True)
        subprocess.run(f"samtools sort -@ {cores} -o {output_bam} {output_file} && samtools index {output_bam}", shell=True, check=True)
        os.remove(output_file)
    return output_bam

def assess_read_coverage(output_bam):
    cov_proportion = {}
    # Run samtools coverage to get coverage statistics
    coverage_cmd = ["samtools", "coverage", output_bam]
    coverage_output = subprocess.check_output(coverage_cmd).decode()
    # Parse coverage output (skip header lines)
    for line in coverage_output.splitlines():
        if line.startswith("#") or line.strip() == "":
            continue
        fields = line.strip().split("\t")
        ref_name = fields[0]
        cov_proportion[ref_name] = float(fields[5])
    return cov_proportion

# make the output directory
if not os.path.exists(args.output):
    os.mkdir(args.output)
if not os.path.exists(os.path.join(args.output, "minimap_out")):
    os.mkdir(os.path.join(args.output, "minimap_out"))
# load the amira TSV
amira_content = pd.read_csv(os.path.join(args.amira_output, "amira_results.tsv"), sep="\t")
# collect fasta sequences
fasta_files = []
for index, row in amira_content.iterrows():
    fasta_file = os.path.join(args.amira_output, "AMR_allele_fastqs", row["Amira allele"], "06.final_sequence.fasta")
    fasta_files.append((row["Amira allele"], fasta_file))
cov_proportion = {}
if not len(fasta_files) == 0:
    for allele, f in tqdm(fasta_files):
        # define the output sam
        output_sam = os.path.join(args.output, "minimap_out", f"{allele}.sam")
        # map the reads to the genes
        output_bam = map_reads_to_fasta(f, args.reads, output_sam, args.cores)
        # get the coverages across the references
        cov_allele = assess_read_coverage(output_bam)
        for a in cov_allele:
            cov_proportion[allele] = cov_allele[a]
# write the results out
outfile = os.path.join(args.output, "amira_coverages.json")
with open(outfile, "w") as o:
    o.write(json.dumps(cov_proportion))
#outfile.to_csv(coverage_df, sep="\t")
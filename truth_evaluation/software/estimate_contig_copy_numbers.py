import subprocess
import pandas as pd
import os
from tqdm import tqdm
import json
import argparse

mapping = {
    "GCA_027944575.1_ASM2794457v1_genomic":  {"ONT": "amira_paper_data/truth_evaluation/nanopore_reads/SRR23044210_1.fastq"},
    "GCA_027944595.1_ASM2794459v1_genomic":  {"ONT": "amira_paper_data/truth_evaluation/nanopore_reads/SRR23044211_1.fastq"},
    "GCA_027944615.1_ASM2794461v1_genomic":  {"ONT": "amira_paper_data/truth_evaluation/nanopore_reads/SRR23044212_1.fastq"},
    "GCA_027944635.1_ASM2794463v1_genomic":  {"ONT": "amira_paper_data/truth_evaluation/nanopore_reads/SRR23044213_1.fastq"},
    "GCA_027944655.1_ASM2794465v1_genomic":  {"ONT": "amira_paper_data/truth_evaluation/nanopore_reads/SRR23044215_1.fastq"},
    "GCA_027944675.1_ASM2794467v1_genomic":  {"ONT": "amira_paper_data/truth_evaluation/nanopore_reads/SRR23044216_1.fastq"},
    "GCA_027944695.1_ASM2794469v1_genomic":  {"ONT": "amira_paper_data/truth_evaluation/nanopore_reads/SRR23044217_1.fastq"},
    "GCA_027944715.1_ASM2794471v1_genomic":  {"ONT": "amira_paper_data/truth_evaluation/nanopore_reads/SRR23044218_1.fastq"},
    "GCA_027944735.1_ASM2794473v1_genomic":  {"ONT": "amira_paper_data/truth_evaluation/nanopore_reads/SRR23044219_1.fastq"},
    "GCA_027944775.1_ASM2794477v1_genomic":  {"ONT": "amira_paper_data/truth_evaluation/nanopore_reads/SRR23044221_1.fastq"},
    "GCA_027944795.1_ASM2794479v1_genomic":  {"ONT": "amira_paper_data/truth_evaluation/nanopore_reads/SRR23044222_1.fastq"},
    "GCA_027944815.1_ASM2794481v1_genomic":  {"ONT": "amira_paper_data/truth_evaluation/nanopore_reads/SRR23044223_1.fastq"},
    "GCA_027944835.1_ASM2794483v1_genomic":  {"ONT": "amira_paper_data/truth_evaluation/nanopore_reads/SRR23044224_1.fastq"},
    "GCA_027944875.1_ASM2794487v1_genomic":  {"ONT": "amira_paper_data/truth_evaluation/nanopore_reads/SRR23044204_1.fastq"},
    "GCA_027944895.1_ASM2794489v1_genomic":  {"ONT": "amira_paper_data/truth_evaluation/nanopore_reads/SRR23044205_1.fastq"},
    "GCA_027944915.1_ASM2794491v1_genomic":  {"ONT": "amira_paper_data/truth_evaluation/nanopore_reads/SRR23044206_1.fastq"},
    "GCA_027944935.1_ASM2794493v1_genomic":  {"ONT": "amira_paper_data/truth_evaluation/nanopore_reads/SRR23044207_1.fastq"},
    "GCA_027944955.1_ASM2794495v1_genomic":  {"ONT": "amira_paper_data/truth_evaluation/nanopore_reads/SRR23044208_1.fastq"},
    "GCA_027945015.1_ASM2794501v1_genomic":  {"ONT": "amira_paper_data/truth_evaluation/nanopore_reads/SRR23044214_1.fastq"},
    "GCA_027945035.1_ASM2794503v1_genomic":  {"ONT": "amira_paper_data/truth_evaluation/nanopore_reads/SRR23044225_1.fastq"},
    "GCA_027945055.1_ASM2794505v1_genomic":  {"ONT": "amira_paper_data/truth_evaluation/nanopore_reads/SRR23044226_1.fastq"},
    "GCA_028551585.1_ASM2855158v1_genomic":  {"ONT": "amira_paper_data/truth_evaluation/nanopore_reads/SRR23044220_1.fastq"},
    "AUSMDU00010405":  {"ONT": "amira_paper_data/truth_evaluation/nanopore_reads/AUSMDU00010405.fastq.gz"},
    "AUSMDU00015264":  {"ONT": "amira_paper_data/truth_evaluation/nanopore_reads/AUSMDU00015264.fastq.gz"},
    "AUSMDU00021208":  {"ONT": "amira_paper_data/truth_evaluation/nanopore_reads/AUSMDU00021208.fastq.gz"},
    "AUSMDU00031899":  {"ONT": "amira_paper_data/truth_evaluation/nanopore_reads/AUSMDU00031899.fastq.gz"},
    "AUSMDU00031978":  {"ONT": "amira_paper_data/truth_evaluation/nanopore_reads/AUSMDU00031978.fastq.gz"},
    "AUSMDU00032793":  {"ONT": "amira_paper_data/truth_evaluation/nanopore_reads/AUSMDU00032793.fastq.gz"},
    "AUSMDU00036400":  {"ONT": "amira_paper_data/truth_evaluation/nanopore_reads/AUSMDU00036400.fastq.gz"},
    "AUSMDU00040126":  {"ONT": "amira_paper_data/truth_evaluation/nanopore_reads/AUSMDU00040126.fastq.gz"},
    "AUSMDU00055259":  {"ONT": "amira_paper_data/truth_evaluation/nanopore_reads/AUSMDU00055259.fastq.gz"},
    "AUSMDU00062512":  {"ONT": "amira_paper_data/truth_evaluation/nanopore_reads/AUSMDU00062512.fastq.gz"}
}

def get_mean_read_depth_per_contig(bam_file):
    # Run samtools depth command and capture output
    result = subprocess.run(
        ["samtools", "depth", bam_file],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=True
    )
    # Read the output into a DataFrame
    data = []
    for line in result.stdout.strip().split("\n"):
        contig, position, depth = line.split("\t")
        data.append((contig, int(position), int(depth)))
    # Create a DataFrame from the parsed data
    df = pd.DataFrame(data, columns=["contig", "position", "depth"])
    # Calculate mean depth for each contig
    mean_depth_per_contig = df.groupby("contig")["depth"].mean().to_dict()
    # Determine the longest contig
    longest_contig = df.groupby("contig")["position"].max().idxmax()
    return mean_depth_per_contig, longest_contig

def map_fastq_to_ref(fastq_file, reference_file, output_sam, cores, minimap2_path="minimap2"):
    # Build the minimap2 command
    output_bam = output_sam.replace(".sam", ".bam")
    command = f"{minimap2_path} --eqx -t {cores} -a -x map-ont --secondary=no -o {output_sam} {reference_file} {fastq_file} && samtools sort -@ {cores} {output_sam} > {output_bam} && rm -rf {output_sam}"
    if not os.path.exists(output_bam):
        subprocess.run(command, shell=True, check=True)
    return output_bam

parser = argparse.ArgumentParser(description="Map genes to assemblies and generate GFF files of alleles with 95% similarity.")
parser.add_argument("assembly", help="Assembly path.")
parser.add_argument("reads", help="Read path.")
parser.add_argument("output", help="Output path.")
parser.add_argument("cores", help="Number of CPUs to use.")
args = parser.parse_args()

sample = os.path.basename(args.assembly).replace(".fasta", "").replace(".fa", "")
sam_file = os.path.join(args.output, sample + ".sam")
reference = args.assembly
output_bam = map_fastq_to_ref(args.reads, args.assembly, sam_file, args.cores)
depth_per_contig, longest_contig = get_mean_read_depth_per_contig(output_bam)
normalised_depths = {c: depth_per_contig[c] / depth_per_contig[longest_contig] for c in depth_per_contig}
with open(os.path.join(args.output, sample + ".json"), "w") as o:
    o.write(json.dumps(normalised_depths))
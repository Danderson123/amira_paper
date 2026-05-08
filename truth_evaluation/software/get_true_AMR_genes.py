import subprocess
import os
import argparse
import pysam

parser = argparse.ArgumentParser()
parser.add_argument("--reads", dest="reads", required=True)
parser.add_argument("--reference", dest="reference", required=True)
parser.add_argument("--output", dest="output", required=True)
parser.add_argument("--cores", dest="cores", required=True)
args = parser.parse_args()

def supplement_with_genes_from_reads(genes_file, nanopore, output_file, cores, similarity_threshold=90.0, length_threshold=100):
    # map the reads to the genes
    cmd = f"minimap2 -a -x map-ont -t {cores} --eqx --MD {genes_file} {nanopore} > {output_file} && samtools sort {output_file} > {output_file.replace('.sam', '.bam')} && samtools index {output_file.replace('.sam', '.bam')}"
    if not os.path.exists(output_file):
        subprocess.run(cmd, shell=True, check=True)
    # import the sam
    valid_candidates = []
    proportion_reference_covered = {}
    max_depth = 0
    with pysam.AlignmentFile(output_file, "r") as sam_file:
        for read in sam_file.fetch():
            if read.is_unmapped:
                continue
            # Calculate matching bases and the length used in the alignment
            matching_bases = 0
            total_aligned_length = 0
            reference_length = sam_file.get_reference_length(read.reference_name)
            for op, length in read.cigartuples:
                if op == 7:  # Match/Mismatch
                    matching_bases += length
                if op != 1 and op != 4 and op != 5:  # Not an insertion to reference
                    total_aligned_length += length
            # Calculate the similarity and the alignment length percentage
            similarity = (matching_bases / total_aligned_length) * 100.0 if total_aligned_length > 0 else 0
            alignment_length_percentage = (total_aligned_length / reference_length) * 100.0 if reference_length > 0 else 0
            # Check against thresholds
            if similarity >= similarity_threshold and alignment_length_percentage >= length_threshold:
                if read.reference_name not in proportion_reference_covered:
                    proportion_reference_covered[read.reference_name] = set()
                proportion_reference_covered[read.reference_name].add(read.query_name)
    for ref in proportion_reference_covered:
        if len(proportion_reference_covered[ref]) > max_depth:
            max_depth = len(proportion_reference_covered[ref])
    for ref in proportion_reference_covered:
        if len(proportion_reference_covered[ref]) >= 5:
            valid_candidates.append((ref, len(proportion_reference_covered[ref])))
    os.remove(output_file)
    return valid_candidates

def apply_rules(gene):
    if "blaCTX-M" in gene:
        gene = "blaCTX-M"
    if "blaNDM" in gene:
        gene = "blaNDM"
    if "blaOXA" in gene:
        gene = "blaOXA"
    if "aac6-Ib" in gene:
        gene = "aac6-Ib"
    if "blaEC" in gene:
        gene = "blaEC"
    # if "oqxA" in gene:
    #     gene = "oqxA"
    # if "oqxB" in gene:
    #     gene = "oqxB"
    if "oqx" in gene or "Oqx" in gene:
        gene = "oqx"
    if "blaTEM" in gene:
        gene = "blaTEM"
    if "fosA" in gene:
        gene = "fosA"
    if "blaCMY" in gene:
        gene = "blaCMY"
    if "aadA" in gene:
        gene = "aadA"
    if "arr-" in gene:
        gene = "arr-"
    if "dfrA" in gene:
        gene = "dfrA"
    if "rmtB" in gene:
        gene = "rmtB"
    if "aac3-II" in gene and "aac3-III" not in gene:
        gene = "aac3-II"
    if "aac3-III" in gene:
        gene = "aac3-III"
    if "blaSHV" in gene:
        gene = "blaSHV"
    if "qnrS" in gene:
        gene = "qnrS"
    if "aac3-I" in gene and "aac3-II" not in gene and "aac3-III" not in gene:
        gene = "aac3-I"
    if "blaKPC" in gene:
        gene = "blaKPC"
    if "mcr-5" in gene:
        gene = "mcr-5"
    if "qnrB" in gene:
        gene = "qnrB"
    if "cmlA" in gene:
        gene = "cmlA"
    if "aph3-I" in gene and "aph3-II" not in gene and "aph3-III" not in gene:
        gene = "aph3-I"
    if "aph3-II" in gene and "aph3-III" not in gene:
        gene = "aph3-II"
    if "aph3-III" in gene:
        gene = "aph(3)-III"
    if "aph3-I" in gene and "aph3-II" not in gene and "aph3-III" not in gene:
        gene = "aph3-I"
    if "aph3-II" in gene and "aph3-III" not in gene:
        gene = "aph3-II"
    if "aph3-III" in gene:
        gene = "aph3-III"
    if "aph4-I" in gene and "aph4-II" not in gene and "aph4-III" not in gene:
        gene = "aph4-I"
    if "aph4-II" in gene and "aph4-III" not in gene:
        gene = "aph4-II"
    if "aph4-III" in gene:
        gene = "aph4-III"
    if "aph6-I" in gene and "aph6-II" not in gene and "aph6-III" not in gene:
        gene = "aph6-I"
    if "aph6-II" in gene and "aph6-III" not in gene:
        gene = "aph6-II"
    if "aph6-III" in gene:
        gene = "aph6-III"
    if "qacE" in gene:
        gene = "qacE"
    if "blaLAP" in gene:
        gene = "blaLAP"
    if "aac6-I" in gene and "aac6-II" not in gene and "aac6-III" not in gene:
        gene = "aac6-I"
    if "aac6-II" in gene and "aac6-III" not in gene:
        gene = "aac6-II"
    if "aac6-III" in gene:
        gene = "aac6-III"
    if "blaDHA" in gene:
        gene = "blaDHA"
    if "qepA" in gene:
        gene = "qepA"
    if "blaIMI" in gene:
        gene = "blaIMI"
    if "mcr-1" in gene:
        gene = "mcr-1"
    return gene

output_sam = os.path.join(args.output, f"{os.path.basename(args.reads).replace('.fastq.gz', '').replace('.gastq', '')}.sam")
hits = supplement_with_genes_from_reads(args.reference, args.reads, output_sam, args.cores)
unique = {}
for h in hits:
    gene, allele = h[0].split(";")
    if gene not in unique:
        unique[gene] = []
    unique[gene].append((gene + ";" + allele, h[1]))
unique_genes = [max(unique[g], key=lambda x: x[1]) for g in unique]
present_genes = [f"{u[0]}\t{u[1]}" for u in sorted(unique_genes, key=lambda x: x[0])]
outfile = os.path.join(args.output, "present_genes.txt")
with open(outfile, "w") as o:
    o.write("\n".join(list(present_genes)))
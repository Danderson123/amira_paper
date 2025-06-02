import subprocess
import os
import argparse
import pysam
import statistics

parser = argparse.ArgumentParser()
parser.add_argument("--reads", dest="reads", required=True)
parser.add_argument("--reference", dest="reference", required=True)
parser.add_argument("--core-genes", dest="core_genes", required=True)
parser.add_argument("--output", dest="output", required=True)
parser.add_argument("--cores", dest="cores", required=True)
args = parser.parse_args()

def supplement_with_genes_from_reads(genes_file, nanopore, output_file, cores, mean_depth, similarity_threshold=90.0, length_threshold=90.0):
    # map the reads to the genes
    cmd = f"minimap2 -a -x map-ont -t {cores} --eqx --MD {genes_file} {nanopore} > {output_file}"
    if not os.path.exists(output_file):
        subprocess.run(cmd, shell=True, check=True)
    # import the sam
    valid_candidates = []
    proportion_reference_covered = {}
    # check if any
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
        #if len(proportion_reference_covered[ref]) >= mean_depth * 0.5:
        valid_candidates.append(ref)
    return valid_candidates

def apply_rules(gene):
    gene = gene.replace("'", "")
    if "blaCTX-M" in gene:
        gene = "blaCTX-M"
    if "blaNDM" in gene:
        gene = "blaNDM"
    if "blaOXA" in gene:
        gene = "blaOXA"
    if "aac(6)-Ib" in gene:
        gene = "aac(6)-Ib"
    if "blaEC" in gene:
        gene = "blaEC"
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
    if "aac(3)-VI" in gene:
        gene = "aac(3)-VI"
    if "aac(3)-II" in gene and "aac(3)-III" not in gene:
        gene = "aac(3)-II"
    if "aac(3)-III" in gene:
        gene = "aac(3)-III"
    if "blaSHV" in gene:
        gene = "blaSHV"
    if "qnrS" in gene:
        gene = "qnrS"
    if "aac(3)-I" in gene and "aac(3)-II" not in gene and "aac(3)-III" not in gene:
        gene = "aac(3)-I"
    if "blaKPC" in gene:
        gene = "blaKPC"
    if "mcr-5" in gene:
        gene = "mcr-5"
    if "mcr-1" in gene:
        gene = "mcr-1"
    if "qnrA" in gene:
        gene = "qnrA"
    if "qnrB" in gene:
        gene = "qnrB"
    if "cmlA" in gene:
        gene = "cmlA"
    if "aph(3)-I" in gene and "aph(3)-II" not in gene and "aph(3)-III" not in gene:
        gene = "aph(3)-I"
    if "aph(3)-II" in gene and "aph(3)-III" not in gene:
        gene = "aph(3)-II"
    if "aph(3)-III" in gene:
        gene = "aph(3)-III"
    if "aph(3)-VI" in gene:
        gene = "aph(3)-VI"
    if "aph(4)-I" in gene and "aph(4)-II" not in gene and "aph(4)-III" not in gene:
        gene = "aph(4)-I"
    if "aph(4)-II" in gene and "aph(4)-III" not in gene:
        gene = "aph(4)-II"
    if "aph(4)-III" in gene:
        gene = "aph(4)-III"
    if "aph(6)-I" in gene and "aph(6)-II" not in gene and "aph(6)-III" not in gene:
        gene = "aph(6)-I"
    if "aph(6)-II" in gene and "aph(6)-III" not in gene:
        gene = "aph(6)-II"
    if "aph(6)-III" in gene:
        gene = "aph(6)-III"
    if "qacE" in gene:
        gene = "qacE"
    if "blaLAP" in gene:
        gene = "blaLAP"
    if "aac(6)-I" in gene and "aac(6)-II" not in gene and "aac(6)-III" not in gene:
        gene = "aac(6)-I"
    if "aac(6)-II" in gene and "aac(6)-III" not in gene:
        gene = "aac(6)-II"
    if "aac(6)-III" in gene:
        gene = "aac(6)-III"
    if "blaDHA" in gene:
        gene = "blaDHA"
    if "qepA" in gene:
        gene = "qepA"
    if "blaIMI" in gene:
        gene = "blaIMI"
    if "ant(2)-I" in gene:
        gene = "ant(2)-I"
    if "ant(3)-I" in gene:
        gene = "ant(3)-I"
    if "ant(9)-I" in gene:
        gene = "ant(9)-I"
    if "blaACC" in gene:
        gene = "blaACC"
    if "blaACT" in gene:
        gene = "blaACT"
    if "blaCARB" in gene:
        gene = "blaCARB"
    if "blaGES" in gene:
        gene = "blaGES"
    if "blaIMP" in gene:
        gene = "blaIMP"
    if "blaGES" in gene:
        gene = "blaGES"
    if "blaVIM" in gene:
        gene = "blaVIM"
    if "toprJ" in gene:
        gene = "toprJ"
    if "blaVEB" in gene:
        gene = "blaVEB"
    if "blaMIR" in gene:
        gene = "blaMIR"
    if "arr" in gene:
        gene = "arr"
    if "mphK" in gene:
        gene = "mph(K)"
    if "satA" in gene:
        gene = "satA"
    if "blaOKP" in gene:
        gene = "blaOKP"
    if "blaHER" in gene:
        gene = "blaHER"
    if "blaMUN" in gene:
        gene = "blaMUN"
    if "blaORN" in gene:
        gene = "blaORN"
    if "catB" in gene:
        gene = "catB"
    if "blaPDC" in gene:
        gene = "blaPDC"
    if "mcr-3" in gene:
        gene = "mcr-3"
    if "catA" in gene:
        gene = "catA"
    if "blaLEN" in gene:
        gene = "blaLEN"
    if "blaPDC" in gene:
        gene = "blaPDC"
    if "blaPDC" in gene:
        gene = "blaPDC"
    if "blaSCO" in gene:
        gene = "blaSCO"
    if "rmtF" in gene:
        gene = "rmtF"
    if "blaFOX" in gene:
        gene = "blaFOX"
    if "blaADC" in gene:
        gene = "blaADC"
    if "blaFRI" in gene:
        gene = "blaFRI"
    if "blaCMH" in gene:
        gene = "blaCMH"
    if "blaSFO" in gene:
        gene = "blaSFO"
    if "cfiA" in gene:
        gene = "cfiA"
    if "ant(6)-I" in gene and not "ant(6)-II" in gene:
        gene = "ant(6)-I"
    if "vanS-B" in gene:
        gene = "vanS-B"
    if "poxtA" in gene:
        gene = "poxtA"
    if "blaR" in gene:
        gene = "blaR"
    if "aadD" in gene:
        gene = "aadD"
    if "mecR" in gene:
        gene = "mecR"
    if "vanC" in gene:
        gene = "vanC"
    return gene

def samtools_get_mean_depth(bam_file):
    # Run samtools coverage and capture output
    result = subprocess.run(
        ["samtools", "coverage", bam_file],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=True,
    )
    # read the output
    mean_depth_per_contig = {}
    for line in result.stdout.strip().split("\n"):
        if line.startswith("#"):
            continue
        rname, startpos, endpos, numreads, covbases, coverage, meandepth, meanbaseq, meanmapqs = (
            line.split("\t")
        )
        if rname not in mean_depth_per_contig:
            mean_depth_per_contig[rname] = 0
        if float(meandepth) > mean_depth_per_contig[rname]:
            mean_depth_per_contig[rname] = float(meandepth)
    return mean_depth_per_contig

def get_core_gene_mean_depth(bam_file, core_gene_file, read_path, cores):
    # map the reads to the core genes
    map_command = f"minimap2 -x map-ont -a --MD --secondary=no -t {cores} "
    map_command += f"-o {bam_file.replace('.bam', '.sam')} "
    map_command += f"{core_gene_file} {read_path} && "
    map_command += (
        f"samtools sort -@ {cores} {bam_file.replace('.bam', '.sam')} > "
    )
    map_command += f"{bam_file} && "
    map_command += f"samtools index {bam_file}"
    subprocess.run(map_command, shell=True, check=True)
    # get the mean depth across each core gene
    mean_depth_per_core_gene = samtools_get_mean_depth(bam_file)
    os.remove(bam_file.replace('.bam', '.sam'))
    os.remove(bam_file)
    os.remove(bam_file + ".bai")
    return statistics.mean(list(mean_depth_per_core_gene.values()))

# define the output sam
output_sam = os.path.join(args.output, f"{os.path.basename(args.reads).replace('.fastq.gz', '').replace('.gastq', '')}.sam")
# get the mean depth across core genes
mean_depth = get_core_gene_mean_depth(output_sam.replace(".sam", ".core_genes.bam"), args.core_genes, args.reads, args.cores)
# call genes found in the reads
hits = supplement_with_genes_from_reads(args.reference, args.reads, output_sam, args.cores, mean_depth)
present_genes = set([apply_rules(h.split(";")[1].split(".NG_")[0]) for h in hits])
outfile = os.path.join(args.output, "present_genes.txt")
with open(outfile, "w") as o:
    o.write("\n".join(list(present_genes)))
import subprocess
import argparse
from Bio import SeqIO
import pyfastaq
import random
import pathlib
import pandas as pd
import os

parser = argparse.ArgumentParser()
parser.add_argument("--context-fasta", dest="context_fasta", required=True)
parser.add_argument("--allele-file", dest="allele_file", required=True)
parser.add_argument("--reference-genome", dest="reference_genome", required=True)
parser.add_argument("--reference-plasmid", dest="reference_plasmid", required=True)
parser.add_argument("--output", dest="output", required=True)
args = parser.parse_args()

# Load the reference assembly
genome = pyfastaq.sequences.file_reader(args.reference_genome)
# Get the genome length and sequence
genome_length, genome_seq = 0, ""
for sequence in genome:
    genome_length = len(sequence.seq)
    genome_seq = str(sequence.seq)

# Load the reference plasmid
plasmid = pyfastaq.sequences.file_reader(args.reference_plasmid)
# Get the plasmid length and sequence
plasmid_length, plasmid_seq = 0, ""
for sequence in plasmid:
    plasmid_length = len(sequence.seq)
    plasmid_seq = str(sequence.seq)

# get the true amr gene content of the block
out_gff = args.output.replace('.fasta', '.AMR_content.gff')
subprocess.run(f"python3 {pathlib.Path(__file__).parent.resolve()}/make_truth_with_minimap.py {args.context_fasta} {args.allele_file} {out_gff}", shell=True, check=True)
# process the gff file
gff_content = {}
with open(out_gff) as i:
    annotations, sequence = i.read().split("##FASTA")
for line in annotations.split("\n"):
    if not line.startswith("#"):
        if line == "":
            continue
        region, source, a_type, start, end, _, strand, dot, feature = line.split("\t")
        if region not in gff_content:
            gff_content[region] = []
        gff_content[region].append((int(start) - 1, int(end) - int(start), feature.replace("Name=", ""), strand))
for r in gff_content:
    gff_content[r] = sorted(gff_content[r], key=lambda x: x[0])
# Load the context file
context = pyfastaq.sequences.file_reader(args.context_fasta)
# Decide where we are putting the AMR blocks
block_insertions = [[], []]
for sequence in context:
    if "_chromosome" in sequence.id:
        block_insertions[0].append((sequence.id, str(sequence.seq)))
    elif "_plasmid" in sequence.id:
        block_insertions[1].append((sequence.id, str(sequence.seq)))
# decide where the blocks will be inserted in the chromosome
chromosome_insertions = sorted(random.sample(range(len(genome_seq)), len(block_insertions[0])))
# decide where the blocks will be inserted in the plasmid
plasmid_insertions = sorted(random.sample(range(len(plasmid_seq)), len(block_insertions[1])))
# decide the plasmid copy number
plasmid_copy_number = random.randint(1, 10)
# split the chromosome into blocks at the insertion points
chromosome_fragments = []
next_start = 0
for position in chromosome_insertions:
    chromosome_fragments.append(genome_seq[next_start:position])
    next_start = position
chromosome_fragments.append(genome_seq[next_start:])
# split the plasmid into blocks at the insertion points
plasmid_fragments = []
next_start = 0
for position in plasmid_insertions:
    plasmid_fragments.append(plasmid_seq[next_start:position])
    next_start = position
plasmid_fragments.append(plasmid_seq[next_start:])
# collect amr gene positions
context_tsv = ["Allele name\tContig\tStart\tEnd\tLength\tStrand\tCopy number"]
# get the chromsome assembly
simulated_chromosome = []
offset = 0
for i in range(len(chromosome_insertions)):
    simulated_chromosome.append(chromosome_fragments[i])
    simulated_chromosome.append(block_insertions[0][i][1])
    adjusted_start = offset + chromosome_insertions[i] + 1
    if block_insertions[0][i][0] in gff_content:
        for amr_start, amr_length, amr_allele, amr_strand in gff_content[block_insertions[0][i][0]]:
            context_tsv.append(f"{amr_allele}\tchromosome\t{adjusted_start + amr_start}\t{adjusted_start + amr_start + amr_length}\t{amr_length}\t{amr_strand}\t1")
    offset += len(block_insertions[0][i][1])
simulated_chromosome.append(chromosome_fragments[-1])
# get the plasmid assembly
simulated_plasmid = []
offset = 0
for i in range(len(plasmid_insertions)):
    simulated_plasmid.append(plasmid_fragments[i])
    simulated_plasmid.append(block_insertions[1][i][1])
    adjusted_start = offset + plasmid_insertions[i] + 1
    if block_insertions[1][i][0] in gff_content:
        for amr_start, amr_length, amr_allele, amr_strand in gff_content[block_insertions[1][i][0]]:
            context_tsv.append(f"{amr_allele}\tplasmid\t{adjusted_start + amr_start}\t{adjusted_start + amr_start + amr_length}\t{amr_length}\t{amr_strand}\t{plasmid_copy_number}")
    offset += len(block_insertions[1][i][1])
simulated_plasmid.append(plasmid_fragments[-1])
# make the simulated assembly
simulated_assembly = []
simulated_chromosome_header = ">contig_1,circular=true"
simulated_assembly.append(simulated_chromosome_header)
simulated_assembly.append("".join(simulated_chromosome))
for i in range(plasmid_copy_number):
    plasmid_header = f">plasmid_{i + 1},circular=true"
    simulated_assembly.append(plasmid_header)
    simulated_assembly.append("".join(simulated_plasmid))
# Join the simulated assembly into a single string
simulated_assembly_str = "\n".join(simulated_assembly)
# write the outputs
if not os.path.exists(os.path.dirname(args.output)):
    os.mkdir(os.path.dirname(args.output))
with open(args.output, "w") as o:
    o.write(simulated_assembly_str)
with open(args.output.replace(".fasta", ".AMR_genes.tsv"), "w") as o:
    o.write("\n".join(context_tsv))
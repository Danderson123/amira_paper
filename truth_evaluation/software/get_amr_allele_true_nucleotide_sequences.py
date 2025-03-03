from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import glob
import io
import os
import json

def parse_gff_fasta(gff_file):
    with open(gff_file) as gff:
        fasta_start = False
        gff_lines = []
        fasta_content = ""
        for line in gff:
            if line.startswith("##FASTA"):
                fasta_start = True
            elif fasta_start:
                fasta_content += line
            else:
                gff_lines.append(line)
    fasta_io = io.StringIO(fasta_content)
    sequences = SeqIO.to_dict(SeqIO.parse(fasta_io, "fasta"))

    gff_df = pd.DataFrame([line.strip().split('\t') for line in gff_lines if not line.startswith('#')],
                          columns=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
    gff_df[['start', 'end']] = gff_df[['start', 'end']].astype(int)
    return gff_df, sequences

def write_gene_sequences(gff_file, output_file, sample_copy_numbers):
    gff_df, sequences = parse_gff_fasta(gff_file)
    with open(output_file, 'w') as output:
        for index, row in gff_df.iterrows():
            if row['type'] == 'CDS':
                #attributes = dict(item.split("=") for item in row['attributes'].split(";"))
                attributes = {"Name": row['attributes'].split("Name=")[1]}
                gene_id = attributes.get('Name', 'Unknown_Name')
                seqid = row['seqid']
                start, end = row['start'], row['end']
                strand = row['strand']
                gene_seq = sequences[seqid].seq[start-1:end]
                if strand == '-':
                    gene_seq = gene_seq.reverse_complement()
                record = SeqRecord(gene_seq, id=f"{gene_id};CCN_{sample_copy_numbers[seqid]}", description="")
                SeqIO.write(record, output, "fasta")

def process_all_gff_files(directory, output_dir, nanopore_contig_copy_numbers):
    for gff_file in glob.glob(directory):
        output_file = os.path.join(output_dir, os.path.basename(gff_file).replace(".gff", ".fasta"))
        sample = os.path.basename(gff_file).replace(".gff", "")
        with open(os.path.join(nanopore_contig_copy_numbers, sample + ".json")) as i:
            sample_copy_numbers = json.load(i)
        write_gene_sequences(gff_file, output_file, sample_copy_numbers)
        print(f'Processed {gff_file} and wrote sequences to {output_file}')

directory = 'truth_jsons/*.gff'  # Replace with your directory containing GFF files
amr_headers = "AMR_alleles_unified.fa"
output_dir = "truth_allele_sequences"
nanopore_contig_copy_numbers = "mapped_nanopore_reads"
with open(amr_headers) as i:
    set_of_amr_genes = set(i.read().split("\n"))
if not os.path.exists(output_dir):
    os.mkdir(output_dir)
process_all_gff_files(directory, output_dir, nanopore_contig_copy_numbers)

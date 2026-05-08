import glob
import os
from tqdm import tqdm

alignment_dir = "/hps/nobackup/iqbal/dander/Amira_panRG_pipeline_test/Escherichia_coli_panRG_thesis/qced_aligned_gene_sequences"
alignments = glob.glob(os.path.join(alignment_dir, "*.fasta"))

all_seqs = []
for f in tqdm(alignments):
	with open(f) as i:
		seqs = i.read().split(">")[1:]
	gene_name = os.path.basename(f).replace(".fasta", "")
	for s in seqs:
		header = s.split("\n")[0]
		cleaned = "".join(s.split("\n")[1:]).replace("-", "")
		all_seqs.append(f">{gene_name};{header}\n{cleaned}")
with open("all_reference_sequences.fasta", "w") as o:
	o.write("\n".join(all_seqs))

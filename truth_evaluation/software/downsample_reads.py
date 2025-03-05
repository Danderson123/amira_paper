import os
import subprocess
import glob
from tqdm import tqdm

all_reads = glob.glob("amira_paper_data/truth_evaluation/nanopore_reads/*")
if not os.path.exists("amira_paper_data/truth_evaluation/subsampled_nanopore_reads"):
    os.mkdir("amira_paper_data/truth_evaluation/subsampled_nanopore_reads")
outdir = "amira_paper_data/truth_evaluation/subsampled_nanopore_reads"
for reads in tqdm(all_reads):
    command = f"rasusa reads --seed 2024 --coverage 200 --genome-size 5mb {reads} > {os.path.join(outdir, os.path.basename(reads).replace('.gz', ''))} && gzip -3 {os.path.join(outdir, os.path.basename(reads).replace('.gz', ''))}"
    subprocess.run(command, shell=True, check=True)
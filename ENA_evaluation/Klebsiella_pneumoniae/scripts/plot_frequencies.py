import glob
import os
from tqdm import tqdm
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist
import sys
from matplotlib.patches import Patch
import json

sys.setrecursionlimit(50000)

amira_outputs = glob.glob("evaluation_results/Amira_output.v0.9.3/*/amira_results.tsv")


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
    if "mcr-8" in gene:
        gene = "mcr-8"
    if "blaCMH" in gene:
        gene = "blaCMH"
    if "blaSFO" in gene:
        gene = "blaSFO"
    if "aacA" in gene:
        gene = "aacA"
    if "cfiA" in gene:
        gene = "cfiA"
    if "qacG" in gene:
        gene = "qacG"
    if "ant(6)-I" in gene and not "ant(6)-II" in gene:
        gene = "ant(6)-I"
    if "aph(7)-I" in gene and "aph(7)-II" not in gene and "aph(7)-III" not in gene:
        gene = "aph(6)-I"
    return gene

with open("AMR_alleles_unified.fa") as i:
    allele_rows = i.read().split(">")[1:]
reference_genes = set()
for r in allele_rows:
    if r != "":
        amira_allele, reference_allele = r.split("\n")[0].split(";")
        reference_genes.add(apply_rules(reference_allele.split(".NG")[0]))

gene_presence_absence = {}
all_samples = set()
gene_in_truth = []
for a in tqdm(amira_outputs):
    sample = os.path.basename(os.path.dirname(a))
    amrfinder = os.path.join("evaluation_results/AMR_finder_plus_results.flye_v2.9.3_nanopore_only_assemblies", sample, "AMR_finder_plus_results.gff")
    verification_file = os.path.join("evaluation_results", "amira_allele_coverages", sample, "amira_coverages.json")
    if not os.path.exists(verification_file):
        continue
    if not os.path.exists(amrfinder):
        continue
    amira_content = pd.read_csv(a, sep="\t")
    amrfinder_content = pd.read_csv(amrfinder, sep="\t")
    with open(verification_file) as i:
        validation_coverages = json.load(i)
    # process amira results
    amira_counts = {}
    for index, row in amira_content.iterrows():
        assert float(validation_coverages[row["Amira allele"]]) >= 85, f"Amira allele {row['Amira allele']} in sample {sample} is below the required coverage threshold ({validation_coverages[row['Amira allele']]})."
        gene_name = apply_rules(row["Determinant name"])
        amira_counts[gene_name] = amira_counts.get(gene_name, 0) + 1
    amrfp_counts = {}
    for _, row in amrfinder_content.iterrows():
        if row["% Coverage of reference sequence"] < 90 or row["Element subtype"] == "POINT":
            continue
        gene_name = apply_rules(row["Gene symbol"])
        if gene_name not in reference_genes:
            continue
        amrfp_counts[gene_name] = amrfp_counts.get(gene_name, 0) + 1

    for gene in set(amira_counts.keys()).union(amrfp_counts.keys()):
        if gene not in gene_presence_absence:
            gene_presence_absence[gene] = {"Amira": 0, "AMRFP Flye": 0}
        if gene in amira_counts:
            gene_presence_absence[gene]["Amira"] += 1
        if gene in amrfp_counts:
            gene_presence_absence[gene]["AMRFP Flye"] += 1
    all_samples.add(sample)

gene_data = []
amira_higher = 0
amira_same = 0
amira_lower = 0

for gene, counts in gene_presence_absence.items():
    amira_freq = counts["Amira"] / len(all_samples)
    amrfp_freq = counts["AMRFP Flye"] / len(all_samples)
    gene_data.append((gene, amira_freq, amrfp_freq))
    if amira_freq > amrfp_freq:
        amira_higher += 1
    elif amira_freq == amrfp_freq:
        amira_same += 1
    else:
        amira_lower += 1

df = pd.DataFrame(gene_data, columns=["Gene", "Amira_Freq", "AMRFP_Flye_Freq"])
df["Max_Freq"] = df[["Amira_Freq", "AMRFP_Flye_Freq"]].max(axis=1)
df = df.sort_values(by="Max_Freq", ascending=False)

plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 16})

fig, ax = plt.subplots(figsize=(36, 12))
ax.set_yscale("log")

def safe_log_value(val, min_val=0.00001):
    return val if val > 0 else min_val

viridis = plt.cm.get_cmap("viridis")
palette = [viridis(i) for i in [0.0, 0.25, 0.5, 0.75, 1.0][::-1]]  # 5 well

for idx, row in df.iterrows():
    x = row["Gene"]
    amira_freq = safe_log_value(row["Amira_Freq"])
    amrfp_freq = safe_log_value(row["AMRFP_Flye_Freq"])

    plt.plot([x, x], [amira_freq, amrfp_freq], color='gray', zorder=1)
    plt.scatter(x, amira_freq, label="Amira", marker="o", color=palette[0], zorder=2,  edgecolor='black', s=200)
    plt.scatter(x, amrfp_freq, label="AMRFP Flye", marker="X",  edgecolor='black', color=palette[1], zorder=2, s=300)
ax.set_xticks(range(len(df)))
ax.set_xticklabels(df["Gene"], rotation=90, fontstyle='italic')
ax.set_ylabel("AMR gene-presence frequency")
ax.set_ylim([0.00001, 1.1])
ax.set_yticks([0.00001, 0.0001, 0.001, 0.01, 0.1, 1])
ax.set_yticklabels(["0", "0.0001", "0.001", "0.01", "0.1", "1"])
ax.set_xlim([-0.5, len(df) - 0.5])

ax.spines['left'].set_visible(True)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_linewidth(0.5)
ax.grid(axis="y", linestyle="--", alpha=0.7, zorder=1)
ax.grid(axis="x", visible=False)
ax.set_axisbelow(True)

plt.tight_layout()
plt.savefig("evaluation_results/gene_frequencies.png", dpi=300)
plt.savefig("evaluation_results/gene_frequencies.pdf")

print(f"Amira frequency higher: {amira_higher} / {len(gene_data)}")
print(f"Amira frequency same: {amira_same} / {len(gene_data)}")
print(f"Amira frequency lower: {amira_lower} / {len(gene_data)}")

print(len(all_samples))
with open("K_pneumoniae_nanopore_accessions.txt", "w") as o:
    o.write("\n".join(list(all_samples)))
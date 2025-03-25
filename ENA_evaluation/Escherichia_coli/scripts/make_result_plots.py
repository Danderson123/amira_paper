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

sys.setrecursionlimit(50000)

amira_outputs = glob.glob("evaluation_results/Amira_v0.8.0_output/*/amira_results.tsv")

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
    return gene

# load all the genes we are interested in
with open("AMR_alleles_unified.fa") as i:
    allele_rows = i.read().split(">")[1:]
reference_genes = set()
ariba_mapping = {}
for r in allele_rows:
    if r != "":
        amira_allele, reference_allele = r.split("\n")[0].split(";")
        reference_genes.add(apply_rules(reference_allele.split(".NG")[0]))

total_counts = {}
pa_data = {"Method": [], "Gene": [], "Val": []}
missing_data = {"sample": [], "gene": [], "component": []}
for a in tqdm(amira_outputs):
    sample = os.path.basename(os.path.dirname(a))
    amrfinder = os.path.join("evaluation_results/AMR_finder_plus_results.flye_v2.9.3_nanopore_only_assemblies", sample, "AMR_finder_plus_results.gff")
    truth_file = os.path.join("evaluation_results", "true_gene_content", sample, "present_genes.txt")
    if not os.path.exists(truth_file):
        continue
    if not os.path.exists(amrfinder):
        continue
    print(sample)
    amira_content = pd.read_csv(a, sep="\t")
    amrfinder_content = pd.read_csv(amrfinder, sep="\t")
    with open(truth_file) as i:
        truth = set([apply_rules(g) for g in i.read().split("\n")])
    # process amira results
    unique_amira_genes = {}
    for index, row in amira_content.iterrows():
        gene_name = apply_rules(row["Determinant name"])
        if gene_name not in unique_amira_genes:
            unique_amira_genes[gene_name] = set()
    # process amrfinder results
    unique_amrfp_genes = set()
    for index, row in amrfinder_content.iterrows():
        if row["% Coverage of reference sequence"] < 90:
            continue
        gene_name = row["Gene symbol"]
        if row["Element subtype"] == "POINT":
            continue
        gene_name = apply_rules(gene_name)
        if gene_name not in reference_genes:
            continue
        unique_amrfp_genes.add(gene_name)
    # detect presence or absence
    for g in truth:
        if g not in total_counts:
            total_counts[g] = set()
        if g in unique_amrfp_genes:
            pa_data["Method"].append("AMRFP Flye")
            pa_data["Gene"].append(g)
            pa_data["Val"].append(1)
        else:
            pa_data["Method"].append("AMRFP Flye")
            pa_data["Gene"].append(g)
            pa_data["Val"].append(0)
        if g in unique_amira_genes:
            pa_data["Method"].append("Amira")
            pa_data["Gene"].append(g)
            pa_data["Val"].append(1)
        else:
            pa_data["Method"].append("Amira")
            pa_data["Gene"].append(g)
            pa_data["Val"].append(0)
        total_counts[g].add(sample)

data = {}
for i in range(len(pa_data["Method"])):
    if pa_data["Method"][i] not in data:
        data[pa_data["Method"][i]] = {}
    if pa_data["Gene"][i] not in data[pa_data["Method"][i]]:
        data[pa_data["Method"][i]][pa_data["Gene"][i]] = 0
    if pa_data["Val"][i] == 1:
        data[pa_data["Method"][i]][pa_data["Gene"][i]] += 1

plot_data = {"Method": [], "Gene": [], "Val": [], "Samples": []}
samples = 0
for m in data:
    for g in data[m]:
        plot_data["Method"].append(m)
        plot_data["Gene"].append(g)
        plot_data["Val"].append(data[m][g] / len(total_counts[g]))
        plot_data["Samples"].append(len(total_counts[g]))
        if len(total_counts[g]) > samples:
            samples = len(total_counts[g])

# Pivot data for heatmap visualization
data_df = pd.DataFrame(plot_data).pivot(index="Method", columns="Gene", values="Val").fillna(0)
# Sort columns by total count (number of samples where the gene is present)
sorted_genes = sorted(data_df.columns, key=lambda gene: (-len(total_counts[gene]), gene))
data_df = data_df[sorted_genes]
data_df = data_df.loc[:, (data_df != 0).any(axis=0)]
# Create custom colormap
cmap = mcolors.LinearSegmentedColormap.from_list("custom_cmap", ["#fde725", "#440154"])

# Plot Normal Heatmap (No clustering)
plt.figure(figsize=(max(8, data_df.shape[1] * 0.55), 3))
ax = sns.heatmap(data_df, cmap=cmap, linewidths=0.5, center=(0 + 1) / 2, vmin=0, vmax=1, 
                 cbar_kws={"label": "Recall"})
plt.xlabel("Gene", fontsize=16, fontname='sans-serif')
plt.ylabel("", fontsize=16, fontname='sans-serif')
ax.set_xticklabels(
    data_df.columns, rotation=90, fontsize=16, fontstyle='italic', fontname='sans-serif'
)
ax.tick_params(axis="y", labelsize=16, rotation=0)
cbar = ax.collections[0].colorbar
cbar.ax.yaxis.label.set_size(16)
cbar.ax.yaxis.label.set_fontname('sans-serif')
for label in cbar.ax.get_yticklabels():
    label.set_fontsize(16)
    label.set_fontname('sans-serif')
plt.tight_layout()
plt.savefig("evaluation_results/amrfp_amira_population_frequencies.pdf")
plt.savefig("evaluation_results/amrfp_amira_population_frequencies.png", dpi=300)
plt.close()
data_df.to_csv("evaluation_results/heatmap_data.tsv", sep="\t")
missing_df = pd.DataFrame(missing_data)
missing_df.to_csv("evaluation_results/missing_genes.tsv", sep="\t", index=False)

# -------------------------------
# Create grouped bar plot for recall
# Each gene will have three bars: one for AMRFP Flye, one for Amira, and one for Amira --no-filtering
scatter_df = pd.DataFrame(plot_data)
# Create a mapping from each Method to a numeric position
methods = scatter_df['Method'].unique()
cat_to_num = {cat: i for i, cat in enumerate(methods)}

# Add a new column with jittered x coordinates:
# Adjust the jitter amplitude (here ±0.1) as needed
scatter_df['x_jitter'] = scatter_df['Method'].map(cat_to_num) + np.random.uniform(-0.25, 0.25, size=len(scatter_df))

plt.figure(figsize=(12.6, 12))
cb_palette = sns.color_palette("viridis", 4)

# Draw the boxplot first (with lower zorder)
ax = sns.boxplot(
    x="Method", 
    y="Val", 
    data=scatter_df, 
    palette=cb_palette[1:],
    whis=1.5,
    fliersize=0,
    zorder=1,
    linewidth=1.5,
    boxprops=dict(alpha=0.9)
)

# Draw the scatter points on top (with a higher zorder)
ax = sns.scatterplot(
    data=scatter_df,
    x="x_jitter",
    y="Val",
    hue="Method",             # color by method
    size="Samples",           # size by sample count
    sizes=(50, 2000),         # Increase marker sizes if needed
    edgecolor="black",
    linewidth=1.5,            # Increase the edge width for better contrast
    palette=cb_palette[1:],
    zorder=2,
    legend=False
)

# Reset x-axis ticks to show the original category labels
ax.set_xticks(list(range(len(methods))))
ax.set_xticklabels(methods, fontsize=16, fontname='sans-serif')

# Compute the minimum, midpoint, and maximum sample counts
min_samples = scatter_df['Samples'].min()
max_samples = scatter_df['Samples'].max()
mid_samples = int((min_samples + max_samples) / 2)

# Define a helper to scale sample counts to marker sizes (using the same mapping as sizes=(50,2000))
def scale_size(s):
    if max_samples == min_samples:
        return 50
    return 50 + (2000 - 50) * ((s - min_samples) / (max_samples - min_samples))

# Create dummy scatter handles for the legend
size_values = [min_samples, max_samples]
size_handles = [
    plt.scatter([], [], s=scale_size(s), color='gray', edgecolor='black', linewidth=1.5)
    for s in size_values
]

# Add the legend (this legend shows only the sample counts)
plt.legend(
    handles=size_handles,
    labels=[str(s) for s in size_values],
    title="Samples",
    labelspacing=1.5,
    fontsize=16,
    title_fontsize=16,
    bbox_to_anchor=(1.01, 0.5),  # adjust these values as needed
    loc="center left",
    borderpad=1,
    frameon=False
)

# Remove y-axis line and vertical grid lines
ax.spines['left'].set_visible(True)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_linewidth(0.5)
ax.grid(axis="y", linestyle="--", alpha=0.7, zorder=1)
ax.grid(axis="x", visible=False)
ax.set_axisbelow(True)

plt.xlabel("", fontsize=16, fontname='sans-serif')
plt.ylabel("Recall", fontsize=16, fontname='sans-serif')
plt.yticks(fontsize=16, fontname='sans-serif')
plt.tight_layout()
# Save the scatter plot
plt.savefig("evaluation_results/amrfp_amira_scatterplot.pdf")
plt.savefig("evaluation_results/amrfp_amira_scatterplot.png", dpi=300)
plt.close()

print("Mean recall for each method:")
print(data_df.mean(axis=1), data_df.std(axis=1))

# Identify genes that occur in only one sample
single_sample_genes = [gene for gene in data_df.columns if len(total_counts[gene]) == 1]

# For each method, compute the proportion of these genes with a recall of 0
for method in data_df.index:
    if single_sample_genes:
        zero_recall = (data_df.loc[method, single_sample_genes] == 1).sum()
        proportion_zero = zero_recall / len(single_sample_genes)
    else:
        proportion_zero = float('nan')
    print(f"Proportion of genes (occurring in 1 sample) with 1 recall for {method}: {proportion_zero:.2f}")

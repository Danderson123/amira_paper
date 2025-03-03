import glob
import os
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import json
import statistics
import pandas as pd
import numpy as np
from matplotlib.ticker import FormatStrFormatter, MultipleLocator
from matplotlib.colors import LinearSegmentedColormap

# Define constants and shared settings
plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 16})

def merge_results(file_list):
    merged_results = {"r9": {}, "r10": {}}
    for f in file_list:
        sample = os.path.basename(f.replace(".json", ""))
        with open(f) as i:
            gene_calls = json.load(i)
        cleaned_calls = {}
        for g in gene_calls:
            gene = g
            gene = apply_rules(gene)
            cleaned_calls[gene] = gene_calls[g]
        if not "AUSM" in sample:
            merged_results["r9"][sample] = cleaned_calls
        else:
            merged_results["r10"][sample] = cleaned_calls
    return merged_results

def process_AMRFP_results(file_list, reference_genes):
    merged_results = {"r9": {}, "r10": {}}
    for f in file_list:
        sample = os.path.basename(os.path.dirname(f))
        if "AUSM" in sample:
            tech = "r10"
        else:
            tech = "r9"
        merged_results[tech][sample] = {}
        results = pd.read_csv(f, sep="\t")
        for index, row in results.iterrows():
            gene = row["Gene symbol"]
            if row["Element subtype"] == "POINT":
                continue
            if "PARTIAL" in row["Method"]:
                continue
            gene = apply_rules(gene)
            if gene not in reference_genes:
                continue
            if gene not in merged_results[tech][sample]:
                merged_results[tech][sample][gene] = 0
            merged_results[tech][sample][gene] += 1
    return merged_results

def apply_rules(gene):
    if "blaCTX-M" in gene:
        gene = "blaCTX-M"
    if "blaNDM" in gene:
        gene = "blaNDM"
    if "blaOXA" in gene:
        gene = "blaOXA"
    if "aac(6')-Ib" in gene:
        gene = "aac(6')-Ib"
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
    if "qnrB" in gene:
        gene = "qnrB"
    if "cmlA" in gene:
        gene = "cmlA"
    if "aph(3'')-I" in gene and "aph(3'')-II" not in gene and "aph(3'')-III" not in gene:
        gene = "aph(3)-I"
    if "aph(3'')-II" in gene and "aph(3'')-III" not in gene:
        gene = "aph(3)-II"
    if "aph(3'')-III" in gene:
        gene = "aph(3)-III"
    if "aph(3')-I" in gene and "aph(3')-II" not in gene and "aph(3')-III" not in gene:
        gene = "aph(3)-I"
    if "aph(3')-II" in gene and "aph(3')-III" not in gene:
        gene = "aph(3)-II"
    if "aph(3')-III" in gene:
        gene = "aph(3)-III"
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
    if "aac(6')-I" in gene and "aac(6')-II" not in gene and "aac(6')-III" not in gene:
        gene = "aac(6')-I"
    if "aac(6')-II" in gene and "aac(6')-III" not in gene:
        gene = "aac(6')-II"
    if "aac(6')-III" in gene:
        gene = "aac(6')-III"
    if "blaDHA" in gene:
        gene = "blaDHA"
    if "qepA" in gene:
        gene = "qepA"
    if "blaIMI" in gene:
        gene = "blaIMI"
    return gene

def process_resfinder_results(resfinder_files, reference_genes):
    import pandas as pd
    resfinder_results = {"r9": {}, "r10": {}}
    for r in resfinder_files:
        sample = os.path.basename(os.path.dirname(r))
        if "AUSM" in sample:
            tech = "r10"
        else:
            tech = "r9"
        results = pd.read_csv(r, sep="\t")
        resfinder_results[tech][sample] = {}
        for index, row in results.iterrows():
            cleaned_gene = apply_rules(row["Resistance gene"])
            if cleaned_gene not in reference_genes:
                continue
            if cleaned_gene not in resfinder_results[tech][sample]:
                resfinder_results[tech][sample][cleaned_gene] = 0
            resfinder_results[tech][sample][cleaned_gene] += 1
    return resfinder_results

def plot_recall_and_precision(truth_results, assembler_results, output):
    # Initialize a list to collect data for plotting
    plot_data = []
    labels = ["Amira", "Flye\nAMRFP", "ResFinder", "Unicycler\nAMRFP"]

    for tech in ["r9", "r10"]:
        for m, method in enumerate(assembler_results):
            label = labels[m]
            recalls = []
            precisions = []
            for sample in truth_results[tech]:
                # Aggregate true positives, false negatives, and false positives across all genes and samples
                total_tp = 0
                total_fn = 0
                total_fp = 0
                if sample not in method[tech]:
                    method[tech][sample] = {}

                for gene in set(truth_results[tech][sample]).union(method[tech].get(sample, {})):
                    truth_count = truth_results[tech][sample].get(gene, 0)
                    method_count = method[tech].get(sample, {}).get(gene, 0)
                    if label == "ResFinder" and method_count > 0:
                        method_count = 1
                    # Calculate true positives, false negatives, and false positives
                    tp = min(truth_count, method_count)
                    fn = max(0, truth_count - method_count)
                    fp = max(0, method_count - truth_count)

                    # Accumulate totals
                    total_tp += tp
                    total_fn += fn
                    total_fp += fp

                # Calculate proportions for stacking
                total_truth_calls = total_tp + total_fn
                total_method_calls = total_tp + total_fp

                tp_truth_proportion = total_tp / total_truth_calls if total_truth_calls > 0 else 0
                fn_truth_proportion = total_fn / total_truth_calls if total_truth_calls > 0 else 0

                tp_method_proportion = total_tp / total_method_calls if total_method_calls > 0 else 0
                fp_method_proportion = total_fp / total_method_calls if total_method_calls > 0 else 0
                recalls.append(tp_truth_proportion)
                precisions.append(tp_method_proportion)

            sensitivity = statistics.mean(recalls)
            fn_prop = 1 - sensitivity
            specificity = statistics.mean(precisions)
            fp_prop = 1- specificity
            print(tech, label, "Recall: ", sensitivity, " Precision: ", specificity, "\n")
            # Append aggregated data for plotting
            plot_data.append({
                "Technology": tech,
                "Method": label,
                "True Positive Proportion (Truth)": sensitivity,
                "False Negative Proportion (Truth)": fn_prop,
                "True Positive Proportion (Method)": specificity,
                "False Positive Proportion (Method)": fp_prop
            })

    # Convert the collected data into a DataFrame for plotting
    df = pd.DataFrame(plot_data)
    plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 16})

    # Create subplots with two rows for r9 and r10
    fig, axes = plt.subplots(2, 2, figsize=(12, 24), sharey=True, sharex=True)

    for row, proportion_type in enumerate([
        ["True Positive Proportion (Truth)", "False Negative Proportion (Truth)"],
        ["True Positive Proportion (Method)", "False Positive Proportion (Method)"]
    ]):
        for col, tech in enumerate(["r9", "r10"]):
            ax = axes[row, col]

            tech_data = df[df["Technology"] == tech]
            methods = tech_data["Method"]
            prop1_values = tech_data[proportion_type[0]]
            prop2_values = tech_data[proportion_type[1]]

            # Plot stacked bars with distinct colors for TP, FN, and FP
            color_tp = "#440d54"  # Dark purple for TP
            color_fn = "#fee724"  # Yellow for FN
            color_fp = "#35b879"  # Green for FP

            if "False Negative" in proportion_type[1]:
                ax.bar(methods, prop1_values, label=proportion_type[0], color=color_tp, zorder=2)
                ax.bar(methods, prop2_values, bottom=prop1_values, label=proportion_type[1], color=color_fn, zorder=2)
            else:
                ax.bar(methods, prop1_values, label=proportion_type[0], color=color_tp, zorder=2)
                ax.bar(methods, prop2_values, bottom=prop1_values, label=proportion_type[1], color=color_fp, zorder=2)

            # Set titles and labels
            row_title = "sensitivity" if row == 0 else "specificity"
            if tech == "r9":
                ax.set_title(f"R9.4.1 {row_title}", fontsize=16)
            if tech == "r10":
                ax.set_title(f"R10.4 {row_title}", fontsize=16)
            if row == 1:
                ax.set_xlabel("Method", fontsize=16)
            if col == 0:
                ax.set_ylabel("Proportion", fontsize=16)
            ax.set_ylim([0.7, 1])
            # Remove y-axis line and vertical grid lines
            ax.spines['left'].set_visible(True)
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.spines['bottom'].set_linewidth(0.5)
            ax.grid(axis="y", linestyle="--", alpha=0.7, zorder=1)
            ax.grid(axis="x", visible=False)
            # Set y-ticks frequency to 0.05
            ax.yaxis.set_major_locator(MultipleLocator(0.05))

    # Add a single legend for the entire figure
    handles = [
        plt.Line2D([0], [0], color=color_tp, lw=10, label="True Positive"),
        plt.Line2D([0], [0], color=color_fn, lw=10, label="False Negative"),
        plt.Line2D([0], [0], color=color_fp, lw=10, label="False Positive")
    ]
    fig.legend(handles=handles, loc="center left", fontsize=16, bbox_to_anchor=(1.005, 0.5), frameon=False)

    # Adjust layout and save the plot
    plt.tight_layout()
    plt.savefig(output, dpi=600)
    plt.savefig(output.replace(".png", ".pdf"))

def plot_cn_heatmap(truth_results, assembler_results, output_prefix):
    import pandas as pd
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    # Initialize a dictionary to hold heatmap data for both technologies
    heatmap_data = {
        "Single": {"Tool": [], "Gene": [], "Recall": []},
        "Multi": {"Tool": [], "Gene": [], "Recall": []}
    }

    tools = ["Amira", "Shovill AMRFP", "Unicycler AMRFP", "Flye AMRFP", "ResFinder"]

    # Iterate through technologies
    for tech in ["r9", "r10"]:
        # Collect per-sample recalls
        for sample, sample_data in truth_results[tech].items():
            for gene, truth_count in sample_data.items():
                if truth_count < 1:
                    continue  # Avoid division by zero or invalid entries

                status = "Single" if truth_count == 1 else "Multi"

                for t, tool_results in enumerate(assembler_results):
                    tool_label = tools[t]
                    tool_count = tool_results.get(tech, {}).get(sample, {}).get(gene, 0)
                    if tool_label == "ResFinder" and tool_count > 0:
                        tool_count = 1
                    recall = tool_count / truth_count  # Calculate recall
                    if recall > 1:
                        recall = 1

                    heatmap_data[status]["Tool"].append(tool_label)
                    heatmap_data[status]["Gene"].append(gene)
                    heatmap_data[status]["Recall"].append(recall)

    # Create DataFrames and calculate mean recall for both single and multi-copy genes
    heatmap_dfs = {}
    pivot_tables = {}
    for status in ["Single", "Multi"]:
        heatmap_dfs[status] = pd.DataFrame(heatmap_data[status])
        pivot_tables[status] = heatmap_dfs[status].groupby(["Tool", "Gene"]).mean().unstack()["Recall"]
    # Define a continuous custom colormap: Yellow (0) -> Purple (1)
    cmap = LinearSegmentedColormap.from_list("RecallMap", ["#fee724", "#440d54"])
    vmin, vmax = 0, 1

    # Generate heatmaps for Single and Multi separately
    for status in ["Single", "Multi"]:
        pivot_table = pivot_tables[status]
        # Adjust figure size dynamically based on the number of columns
        fig_width = max(8, pivot_table.shape[1] * 0.55)  # Adjust scaling factor as needed
        fig, ax = plt.subplots(figsize=(fig_width, 5))

        sns.heatmap(
            pivot_table,
            cmap=cmap,
            center=(vmin + vmax) / 2,
            annot=False,
            cbar_kws={"label": "Mean Recall"},
            ax=ax,
            vmin=vmin,
            vmax=vmax
        )
        ax.set_title(f"{status}-copy gene mean recall", fontsize=16)
        ax.set_xlabel("Gene", fontsize=16)
        ax.set_ylabel("", fontsize=16)
        ax.set_xticklabels(
            pivot_table.columns, rotation=90, fontsize=16, fontstyle='italic'
        )
        ax.tick_params(axis="y", labelsize=16, rotation=0)

        # Adjust layout
        plt.tight_layout()

        # Save the heatmap
        output_file = f"{output_prefix}_{status.lower()}_heatmap.png"
        plt.savefig(output_file, dpi=600)
        plt.savefig(output_file.replace(".png", ".pdf"))
        plt.close()


truth_dir = "truth_jsons"
amira_dir = "amira_output"
flye_dir = "AMR_finder_plus_results.flye_v2.9.3_nanopore_only_assemblies"
unicycler_dir = "AMR_finder_plus_results.unicycler_v0.5.0_hybrid_assemblies"
resfinder_dir = "resfinder_results"
shovill_dir = "AMR_finder_plus_results.shovill_v1.1.0_illumina_only_assemblies"
allele_file = "AMR_alleles_unified.fa"
output_dir = "truth_results"

if not os.path.exists(output_dir):
    os.mkdir(output_dir)
# load all the genes we are interested in
with open(allele_file) as i:
    allele_rows = i.read().split(">")[1:]
reference_genes = set()
for r in allele_rows:
    if r != "":
        amira_allele, reference_allele = r.split("\n")[0].split(";")
        reference_genes.add(apply_rules(reference_allele.split(".NG")[0]))
# merge the results
truth_results = merge_results(glob.glob(os.path.join(truth_dir, "*.json")))
flye_results = process_AMRFP_results(glob.glob(os.path.join(flye_dir, "*", "*.tsv")), reference_genes)
unicycler_results = process_AMRFP_results(glob.glob(os.path.join(unicycler_dir, "*", "*.tsv")), reference_genes)
resfinder_results = process_resfinder_results(glob.glob(os.path.join(resfinder_dir, "*", "ResFinder_results_tab.txt")), reference_genes)
shovill_results = process_AMRFP_results(glob.glob(os.path.join(shovill_dir, "*", "*.tsv")), reference_genes)
# process the amira results
amira_results = {"r9": {}, "r10": {}}
samples = []
for s in glob.glob(os.path.join(amira_dir, "*")):
    if os.path.exists(os.path.join(s, "amira_results.tsv")):
        sample = os.path.basename(s)
        samples.append(sample)
        amira_table = pd.read_csv(os.path.join(s, "amira_results.tsv"), sep="\t")
        if "AUSM" in sample:
            amira_results["r10"][sample] = {}
            r10 = True
        else:
            amira_results["r9"][sample] = {}
            r10 = False
        for index, row in amira_table.iterrows():
            reference_gene = apply_rules(row["Determinant name"])
            if r10 is True:
                if reference_gene not in amira_results["r10"][sample]:
                    amira_results["r10"][sample][reference_gene] = 0
                amira_results["r10"][sample][reference_gene] += 1
            if r10 is False:
                if reference_gene not in amira_results["r9"][sample]:
                    amira_results["r9"][sample][reference_gene] = 0
                amira_results["r9"][sample][reference_gene] += 1
# compensate for structural variants that we are going to ignore
if "AUSMDU00021208" in amira_results["r10"]:
    amira_results["r10"]["AUSMDU00021208"][apply_rules("blaTEM-1")] = amira_results["r10"]["AUSMDU00021208"][apply_rules("blaTEM-1")] - 1
if "GCA_027944615.1_ASM2794461v1_genomic" in amira_results["r9"]:
    amira_results["r9"]["GCA_027944615.1_ASM2794461v1_genomic"][apply_rules("blaTEM-116")] = amira_results["r9"]["GCA_027944615.1_ASM2794461v1_genomic"][apply_rules("blaTEM-116")] - 1

# plot the recall and precisions of each tool
plot_recall_and_precision(truth_results,
                    [
                        amira_results,
                        flye_results,
                        resfinder_results,
                        unicycler_results
                    ],
                    os.path.join(output_dir, "figure_4a.png"))
# plot the recall of each single and multi copy AMR gene
plot_cn_heatmap(truth_results,
                    [
                        amira_results,
                        shovill_results,
                        unicycler_results,
                        flye_results,
                        resfinder_results,
                    ],
                    os.path.join(output_dir, "figure_4d_cn_heatmap.png"))

# Initialize dictionary to store data per scenario, depth, and length
data_list = {s: {} for s in samples}

def rc_sequence(sequence):
    replacement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(list(reversed([replacement[b] for b in list(sequence)])))

all_similarities = {}
all_copy_number_tuples = {}
for s in tqdm(samples):
    if "AUSM" in s:
        method = "10.4"
    else:
        method = "9.4.1"
    if method not in all_similarities:
        all_similarities[method] = []
        all_copy_number_tuples[method] = []
    # Import the truth json
    with open(os.path.join("truth_jsons", s + ".json")) as i:
        truth_counts = json.load(i)
    # import the truth fasta
    with open(os.path.join("truth_allele_sequences", f"{s}.fasta")) as i:
        allele_seqs = i.read().split(">")[1:]
    true_nucleotide_sequences = {}
    true_copy_numbers = {}
    for allele in allele_seqs:
        allele_name = allele.split("\n")[0]
        truth_sequence = "".join(allele.split("\n")[1:])
        amira_allele, reference_allele, cellular_copy_number = allele_name.split(";")
        gene_name = apply_rules(reference_allele.split(".")[0])
        if gene_name not in true_nucleotide_sequences:
            true_nucleotide_sequences[gene_name] = []
            true_copy_numbers[gene_name] = []
        true_nucleotide_sequences[gene_name].append(f">{gene_name}_truth\n{truth_sequence}")
        true_copy_numbers[gene_name].append(float(cellular_copy_number.replace("CCN_", "")))
    # load the amira tsv
    amira_out = os.path.join("amira_output", s, "amira_results.tsv")
    if not os.path.exists(amira_out):
        continue
    amira_results = pd.read_csv(amira_out, sep="\t")
    # get the amira allele sequences
    amira_nucleotide_sequences = {}
    amira_copy_numbers = {}
    for index, row in amira_results.iterrows():
        amira_allele_name = row["Amira allele"]
        gene_name = apply_rules(row["Determinant name"])
        with open(os.path.join("amira_output", s, "AMR_allele_fastqs", amira_allele_name, "06.final_sequence.fasta")) as i:
            allele_sequence = "".join(i.read().split("\n")[1:])
        if gene_name not in amira_nucleotide_sequences:
            amira_nucleotide_sequences[gene_name] = []
            amira_copy_numbers[gene_name] = []
        amira_nucleotide_sequences[gene_name].append(f">{gene_name}_amira\n{allele_sequence}")
        amira_copy_numbers[gene_name].append(row['Approximate copy number'])
    # get the nucleotide accuracy of the amira alleles
    for gene_name in true_nucleotide_sequences:
        if gene_name in amira_nucleotide_sequences:
            all_sequences = "\n".join(true_nucleotide_sequences[gene_name] + amira_nucleotide_sequences[gene_name])
            #amira_similarities = calculate_allele_accuracy_with_mafft(all_sequences, output_dir)
            amira_similarities, copy_number_tuples = calculate_allele_accuracy_with_mafft(all_sequences, output_dir, true_copy_numbers[gene_name], amira_copy_numbers[gene_name])
            all_similarities[method] += amira_similarities
            all_copy_number_tuples[method] += copy_number_tuples

#plot_nucleotide_results(all_similarities, os.path.join(output_dir, "nucleotide_accuracies.png"))
plot_nucleotide_results_violin(all_similarities, os.path.join(output_dir, "figure_4c.png"))
plot_copy_numbers(all_copy_number_tuples, os.path.join(output_dir, "figure_4b.png"))
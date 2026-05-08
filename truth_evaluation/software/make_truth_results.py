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
import subprocess
from Bio import AlignIO
from matplotlib.lines import Line2D

from get_assembler_allele_accuracies import calculate_assembler_accuracies

# Define constants and shared settings
plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 16})

def merge_results(file_list):
    merged_results = {"r9": {}, "r10": {}}
    truth_alleles =  {"r9": {}, "r10": {}}
    for f in file_list:
        sample = os.path.basename(f.replace(".json", ""))
        with open(f) as i:
            gene_calls = json.load(i)
        cleaned_calls = {}
        sample_alleles = {}
        for g in gene_calls:
            gene = g
            gene = apply_rules(gene)
            cleaned_calls[gene] = gene_calls[g]
        with open(f.replace(".json", ".allele_names.json")) as i:
            allele_names = json.load(i)
        for g in allele_names:
            if g not in sample_alleles:
                sample_alleles[g] = 0
            sample_alleles[g] += 1
        if not "AUSM" in sample:
            merged_results["r9"][sample] = cleaned_calls
            truth_alleles["r9"][sample] = sample_alleles
        else:
            merged_results["r10"][sample] = cleaned_calls
            truth_alleles["r10"][sample] = sample_alleles
    return merged_results, truth_alleles

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
    if "mcr-1" in gene:
        gene = "mcr-1"
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
            if row["Identity"] < 90 or row["Coverage"] < 90:
                continue
            if cleaned_gene not in resfinder_results[tech][sample]:
                resfinder_results[tech][sample][cleaned_gene] = 0
            resfinder_results[tech][sample][cleaned_gene] += 1
    return resfinder_results

def plot_recall_and_precision(truth_results, assembler_results, output):
    # Initialize a list to collect data for plotting
    plot_data = []
    labels = ["Amira", "Flye AMRFP", "Raven AMRFP", "Wtdbg2 AMRFP", "ResFinder", "Autocycler AMRFP"]
    total = {l: [] for l in labels}
    for tech in ["r9", "r10"]:
        for m, method in enumerate(assembler_results):
            label = labels[m]
            recalls = []
            precisions = []
            counts = {}
            for sample in truth_results[tech]:
                # Aggregate true positives, false negatives, and false positives across all genes and samples
                total_tp = 0
                total_fn = 0
                total_fp = 0
                if sample not in method[tech]:
                    continue
                    method[tech][sample] = {}
                counts[sample] = {}
                for gene in set(truth_results[tech][sample]).union(method[tech].get(sample, {})):
                    truth_count = truth_results[tech][sample].get(gene, 0)
                    method_count = method[tech].get(sample, {}).get(gene, 0)
                    counts[sample][gene] = {"truth": truth_count, label: method_count}
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
                if "Amira" in label and total_fp > 0:
                    print(sample, tech, truth_results[tech][sample], method[tech][sample], [k for k in method[tech][sample] if method[tech][sample].get(k, 0) < truth_results[tech][sample].get(k, 0)])
                precisions.append(tp_method_proportion)
            total[label] += precisions
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
    plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 20})

    # Create one plot per technology (R9 and R10)
    fig, ax = plt.subplots(figsize=(12, 12))

    # Sample distinct points across the viridis colormap
    viridis = plt.cm.get_cmap("viridis")
    palette = [viridis(i) for i in [0.0, 0.166, 0.333, 0.5, 0.666, 0.833, 1.0][::-1]]  # 7 well-separated colors
    methods = df["Method"].unique()
    techs = ["R9.4.1", "R10.4.1"]

    # Assign color to method, shape to technology
    marker_shapes = {"R9.4.1": "d", "R10.4.1": "o"}
    method_colors = {method: palette[i % len(palette)] for i, method in enumerate(methods)}

    for _, row in df.iterrows():
        tech = row["Technology"]
        if tech == "r9":
            tech = "R9.4.1"
        if tech == "r10":
            tech = "R10.4.1"
        method = row["Method"]
        ax.scatter(
            row["True Positive Proportion (Truth)"], row["True Positive Proportion (Method)"],
            s=500,
            color=method_colors[method],
            marker=marker_shapes[tech],
            edgecolor='black',
            label=f"{method}_{tech}",
            zorder=0
        )

    # Create legend for methods (colors)
    method_handles = [
        Line2D([0], [0], marker='s', color='w', label=method,
            markerfacecolor=color, markeredgecolor="black", markersize=15, linestyle='None')
        for method, color in method_colors.items()
    ]

    # Create legend for technologies (shapes)
    tech_handles = [
        Line2D([0], [0], marker=marker_shapes[tech], color='black', label=f"R{tech[1:]}", 
            markerfacecolor='grey', markersize=15, linestyle='None')
        for tech in techs
    ]

    # Combine and place legend at bottom right
    all_handles = tech_handles + method_handles
    ax.legend(handles=all_handles, loc='lower right', fontsize=20, title_fontsize=20, frameon=True)
     # Configure plot
    ax.set_xlabel("Genomic-copy-number recall", fontsize=20)
    ax.set_ylabel("Genomic-copy-number precision", fontsize=20)
    ax.set_xlim([0.69, 1.01])
    ax.set_ylim([0.69, 1.01])
    ax.xaxis.set_major_locator(MultipleLocator(0.05))
    ax.yaxis.set_major_locator(MultipleLocator(0.05))
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.gca().spines['bottom'].set_linewidth(1)
    plt.gca().spines['left'].set_linewidth(1)
    plt.gca().spines['bottom'].set_color("black")
    plt.gca().spines['left'].set_color("black")
    ax.grid(axis="y", linestyle="--", alpha=0.7, zorder=1)
    ax.grid(axis="x", linestyle="--", alpha=0.7, zorder=1)

    plt.tight_layout()
    plt.savefig(output, dpi=600)
    plt.savefig(output.replace(".png", ".pdf"))

def print_single_and_multicopy_stats(truth_results, assembler_results, tools):
    print("\n===== Mean Recall and Precision by Tech, Tool, and Copy Type =====")
    all_precision_data = {status: {tool: [] for tool in tools} for status in ["Single", "Multi"]}
    all_recall_data = {status: {tool: [] for tool in tools} for status in ["Single", "Multi"]}
    # Also store overall (across Single+Multi) aggregated across tech
    overall_precision_data = {tool: [] for tool in tools}
    overall_recall_data = {tool: [] for tool in tools}
    for tech in ["r9", "r10"]:
        precision_data = {status: {tool: [] for tool in tools} for status in ["Single", "Multi"]}
        recall_data = {status: {tool: [] for tool in tools} for status in ["Single", "Multi"]}
        for sample, sample_data in truth_results[tech].items():
            for gene, truth_count in sample_data.items():
                status = "Single" if truth_count == 1 else "Multi"
                for t, tool_results in enumerate(assembler_results):
                    tool_label = tools[t]
                    method_count = tool_results.get(tech, {}).get(sample, {}).get(gene, 0)
                    tp = min(truth_count, method_count)
                    fn = max(0, truth_count - method_count)
                    fp = max(0, method_count - truth_count)

                    total_truth_calls = tp + fn
                    total_method_calls = tp + fp

                    tp_truth_proportion = tp / truth_count if truth_count > 0 else 0
                    tp_method_proportion = tp / method_count if method_count > 0 else 0

                    recall_data[status][tool_label].append(tp_truth_proportion)
                    precision_data[status][tool_label].append(tp_method_proportion)

                    all_recall_data[status][tool_label].append(tp_truth_proportion)
                    all_precision_data[status][tool_label].append(tp_method_proportion)

                    overall_recall_data[tool_label].append(tp_truth_proportion)
                    overall_precision_data[tool_label].append(tp_method_proportion)
    for status in ["Single", "Multi"]:
            print(f"\nTech: {tech.upper()}, Gene Type: {status}")
            for tool in tools:
                mean_recall = np.mean(recall_data[status][tool]) * 100 if recall_data[status][tool] else 0
                mean_precision = np.mean(precision_data[status][tool]) * 100 if precision_data[status][tool] else 0
                print(f"  {tool:18s} | Recall: {mean_recall:6.2f}% | Precision: {mean_precision:6.2f}%")
    for status in ["Single", "Multi"]:
        print(f"\nTech: Aggregated, Gene Type: {status}")
        for tool in tools:
            mean_recall = np.mean(all_recall_data[status][tool]) * 100 if all_recall_data[status][tool] else 0
            mean_precision = np.mean(all_precision_data[status][tool]) * 100 if all_precision_data[status][tool] else 0
            print(f"  {tool:18s} | Recall: {mean_recall:6.2f}% | Precision: {mean_precision:6.2f}%")
    print("\n===== Overall (All Techs, All Copy Types) =====")
    for tool in tools:
        mean_recall = np.mean(overall_recall_data[tool]) * 100 if overall_recall_data[tool] else 0
        mean_precision = np.mean(overall_precision_data[tool]) * 100 if overall_precision_data[tool] else 0
        print(f"  {tool:18s} | Recall: {mean_recall:6.2f}% | Precision: {mean_precision:6.2f}%")

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

    tools = ["Amira", "Flye AMRFP", "Raven AMRFP", "Wtdbg2 AMRFP", "ResFinder", "Autocycler AMRFP"]

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
                    # if tool_label == "ResFinder" and tool_count > 1:
                    #     print(sample, gene)
                    #     tool_count = 1
                    recall = tool_count / truth_count  # Calculate recall
                    if recall > 1:
                        recall = 1

                    heatmap_data[status]["Tool"].append(tool_label)
                    heatmap_data[status]["Gene"].append(gene)
                    heatmap_data[status]["Recall"].append(recall)
    # print the mean single and multicopy gene stats
    print_single_and_multicopy_stats(truth_results, assembler_results, tools)
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
        # Reorder rows based on external list
        pivot_table = pivot_table.reindex(tools)
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
        ax.set_title(f"{status}-copy gene mean recall", fontsize=20)
        ax.set_xlabel("Gene", fontsize=20)
        ax.set_ylabel("", fontsize=20)
        ax.set_xticklabels(
            pivot_table.columns, rotation=90, fontsize=20, fontstyle='italic'
        )
        ax.tick_params(axis="y", labelsize=20, rotation=0)

        # Adjust layout
        plt.tight_layout()

        # Save the heatmap
        output_file = f"{output_prefix}_{status.lower()}_heatmap.png"
        plt.savefig(output_file, dpi=600)
        plt.savefig(output_file.replace(".png", ".pdf"))
        plt.close()

def calculate_allele_accuracy_with_mafft(all_seqs, output_dir, true_c_n, amira_c_n):
    if not os.path.exists(os.path.join(output_dir, "temp_files")):
        os.mkdir(os.path.join(output_dir, "temp_files"))
    # Create a combined fasta file
    combined_fasta = os.path.join(output_dir, "temp_files", "combined.fasta")
    with open(combined_fasta, "w") as combined:
        combined.write(all_seqs)
    # Run MAFFT on the combined fasta file
    mafft_command = ["mafft", "--auto", "--quiet", combined_fasta]
    aligned_fasta = combined_fasta.replace(".fasta", ".aligned.fasta")
    with open(aligned_fasta, "w") as aligned:
        subprocess.run(mafft_command, stdout=aligned)
    # Load the alignment
    alignment = AlignIO.read(aligned_fasta, "fasta")
    # Extract sequences
    seqs = [(record.id, str(record.seq).upper()) for record in alignment]
    truth_seqs = [aligned for header, aligned in seqs if "_truth" in header]
    amira_seqs = [aligned for header, aligned in seqs if "_amira" in header]
    # Create a similarity matrix
    similarity_matrix = np.zeros((len(truth_seqs), len(amira_seqs)))
    # Fill the similarity matrix
    for i, truth_seq in enumerate(truth_seqs):
        for j, amira_seq in enumerate(amira_seqs):
            matching = 0
            gapless = 0
            for b in range(len(truth_seq)):
                #if truth_seq[b] != "-" and amira_seq[b] != "-":
                if truth_seq[b] == amira_seq[b]:
                    matching += 1
                gapless += 1
            similarity = matching / gapless if gapless > 0 else 0
            similarity_matrix[i, j] = similarity
    # Perform the pairing
    paired_similarities = []
    paired_truths = set()
    paired_amiras = set()
    cn_tuples = []
    while len(paired_truths) < len(truth_seqs) and len(paired_amiras) < len(amira_seqs):
        # Find the highest similarity in the matrix that hasn't been paired yet
        max_similarity = -1
        best_truth_idx = -1
        best_amira_idx = -1
        copy_number_similarity = 100000
        for i in range(len(truth_seqs)):
            if i in paired_truths:
                continue
            for j in range(len(amira_seqs)):
                if j in paired_amiras:
                    continue
                if similarity_matrix[i, j] > max_similarity:
                #if abs(true_c_n[i] - amira_c_n[j]) <= copy_number_similarity:
                    max_similarity = similarity_matrix[i, j]
                    best_truth_idx = i
                    best_amira_idx = j
                    copy_number_similarity = abs(true_c_n[i] - amira_c_n[j])
                if similarity_matrix[i, j] == max_similarity:
                    copy_number_diff = abs(true_c_n[i] - amira_c_n[j])
                    if copy_number_diff < copy_number_similarity:
                        max_similarity = similarity_matrix[i, j]
                        best_truth_idx = i
                        best_amira_idx = j
                        copy_number_similarity = abs(true_c_n[i] - amira_c_n[j])
        # If a valid pair was found, mark the truth and amira alleles as paired
        if best_truth_idx != -1 and best_amira_idx != -1:
            paired_similarities.append(max_similarity)
            cn_tuples.append((true_c_n[best_truth_idx], amira_c_n[best_amira_idx]))
            paired_truths.add(best_truth_idx)
            paired_amiras.add(best_amira_idx)
    return paired_similarities, cn_tuples

def old_plot_nucleotide_results_violin(similarities, output_file):
    combined_data = []
    all_accuracies = []

    # Collect data + print means
    for method_name, tech_dict in similarities.items():
        for tech in ["9.4.1", "10.4.1"]:
            data = tech_dict.get(tech, [])
            tech_label = "R" + tech
            for value in data:
                combined_data.append({
                    "Method": method_name,
                    "Technology": tech_label,
                    "Allele Accuracy": float(value)
                })
                all_accuracies.append(float(value))
            if data:
                print(f"{method_name} - {tech_label} accuracy = {statistics.mean(data):.4f}")

    if all_accuracies:
        print(f"All allele accuracies = {statistics.mean(all_accuracies):.4f}")
    else:
        print("No allele accuracies found.")
        return

    df = pd.DataFrame(combined_data)

    # Order methods by overall mean (cleaner to read)
    method_order = (
        df.groupby("Method")["Allele Accuracy"]
          .mean()
          .sort_values(ascending=False)
          .index
          .tolist()
    )
    df["Method Label"] = pd.Categorical(df["Method"], categories=method_order, ordered=True)

    # Consistent hue order
    hue_order = ["R9.4.1", "R10.4.1"]
    present_hues = [h for h in hue_order if h in set(df["Technology"])]

    plt.figure(figsize=(12, 12), dpi=600)
    ax = plt.gca()

    sns.violinplot(
        y="Method Label",
        x="Allele Accuracy",
        hue="Technology",
        hue_order=present_hues,
        data=df,
        inner="quartile",
        linewidth=1.2,
        palette="colorblind",
        cut=0,
        bw_adjust=0.9,
        dodge=True,
        ax=ax
    )

    sns.stripplot(
        y="Method Label",
        x="Allele Accuracy",
        hue="Technology",
        hue_order=present_hues,
        data=df,
        dodge=True,
        jitter=0.18,
        marker="o",
        edgecolor="black",
        linewidth=0.6,
        alpha=0.65,
        size=5.5,
        palette="colorblind",
        legend=False,
        ax=ax
    )

    # --- Styling (0 to 1 as requested) ---
    ax.set_xlim(0.0, 1.0)

    # Reference lines
    ax.axvline(1.0, linestyle="--", linewidth=1.0, alpha=0.6)
    for v in [0.90, 0.95, 0.98, 0.99]:
        ax.axvline(v, linestyle=":", linewidth=0.8, alpha=0.30)

    # Ticks: coarse across 0-1 + extra near 1
    ticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 0.90, 0.95, 0.98, 0.99]
    ticks = sorted(set(ticks))
    ax.set_xticks(ticks)

    ax.set_xlabel("Allele Accuracy", fontsize=18, labelpad=10)
    ax.set_ylabel("")
    ax.tick_params(axis="x", labelsize=14)
    ax.tick_params(axis="y", labelsize=14)

    # Clean spines + subtle grids
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_linewidth(1.0)
    ax.spines["bottom"].set_linewidth(1.0)

    ax.grid(axis="x", linestyle="--", alpha=0.30)
    ax.grid(axis="y", linestyle=":", alpha=0.20)

    # Legend: dedupe entries (violin+strip causes repeats)
    handles, labels = ax.get_legend_handles_labels()
    if len(labels) >= len(present_hues):
        handles = handles[:len(present_hues)]
        labels = labels[:len(present_hues)]
    ax.legend(
        handles,
        labels,
        loc="lower right",
        fontsize=14,
        title="Technology",
        title_fontsize=14,
        frameon=True
    )

    plt.tight_layout()
    plt.savefig(output_file, dpi=900)
    if output_file.lower().endswith(".png"):
        plt.savefig(output_file[:-4] + ".pdf")
    else:
        plt.savefig(output_file + ".pdf")
    plt.close()

def plot_nucleotide_results_violin(similarities, output_file):
    combined_data = []
    all_accuracies = []

    for method_name, tech_dict in similarities.items():
        for tech in ["9.4.1", "10.4.1"]:
            data = tech_dict.get(tech, [])
            tech = "R" + tech
            for value in data:
                combined_data.append({
                    "Method": method_name,
                    "Technology": tech,
                    "Allele Accuracy": value
                })
                all_accuracies.append(value)
            if data:
                print(f"{method_name} - {tech} accuracy = {statistics.mean(data):.4f}")

    print(f"All allele accuracies = {statistics.mean(all_accuracies):.4f}")

    # Convert to DataFrame
    df = pd.DataFrame(combined_data)

    # Create a combined label for x-axis (one per method)
    df["Method Label"] = df["Method"]

    # Set up the figure
    plt.figure(figsize=(24, 12), dpi=600)

    # Split violin plot: plot both technologies side-by-side within each method
    sns.violinplot(
        x="Method Label",
        y="Allele Accuracy",
        hue="Technology",
        data=df,
        split=False,
        inner="quartile",
        linewidth=1.5,
        palette="colorblind",
        cut=0
    )

    # Add stripplot over violins for individual points
    sns.stripplot(
        x="Method Label",
        y="Allele Accuracy",
        hue="Technology",
        data=df,
        dodge=True,
        jitter=True,
        marker='o',
        edgecolor='black',
        linewidth=1,
        alpha=1,
        size=10,
        palette="colorblind",
        legend=False
    )
    # Styling
    plt.ylabel("Allele Accuracy", fontsize=20, labelpad=15)
    plt.xlabel("", fontsize=20)
    plt.ylim([0.6, 1])
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    ax = plt.gca()
    ax.spines['left'].set_visible(True)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.gca().spines['bottom'].set_linewidth(1)
    plt.gca().spines['left'].set_linewidth(1)
    plt.gca().spines['bottom'].set_color("black")
    plt.gca().spines['left'].set_color("black")
    ax.grid(axis="y", linestyle="--", alpha=0.7, zorder=1)
    #ax.grid(axis="x", visible=False)
    ax.legend(loc='center right', fontsize=20, title_fontsize=20, frameon=True)
    plt.tight_layout()
    plt.savefig(output_file, dpi=900)
    plt.savefig(output_file.replace(".png", ".pdf"))
    plt.close()

def plot_copy_numbers(copy_number_tuples_by_depth, output_file):
    plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 20})
    fig = plt.figure(figsize=(12, 12), dpi=600)  # Explicitly set DPI here
    line_styles = ['-', '--', '-.', ':']  # Line styles for the lines
    marker_styles = ['o', 's', '^', 'D']  # Markers for the points
    x_vals = []
    y_vals = []
    for i, d in enumerate(copy_number_tuples_by_depth):
        x_vals += [float(t[0]) for t in copy_number_tuples_by_depth[d]]
        y_vals += [float(t[1]) for t in copy_number_tuples_by_depth[d]]
    # Compute R² with respect to the identity line y = x
    x_vals = np.array(x_vals)
    y_vals = np.array(y_vals)
    y_mean = np.mean(y_vals)
    ss_tot = np.sum((y_vals - y_mean) ** 2)
    ss_res = np.sum((y_vals - x_vals) ** 2)
    r2 = 1 - (ss_res / ss_tot)
    print(f"R² value (vs identity line y = x): {r2:.4f}")
    # Sample distinct points across the viridis colormap
    viridis = plt.cm.get_cmap("viridis")
    palette = [viridis(i) for i in [0.0, 0.25, 0.5, 0.75, 1.0][::-1]]  # 5 well
    # Scatter plot for data points
    plt.scatter(x_vals, y_vals, s=60, color=palette[0], edgecolor='black')
    # Plot a reference line
    plt.plot([i for i in range(9)], [i for i in range(9)], linestyle="--", color="darkgrey")
    # Styling the plot
    plt.grid(axis="y", linestyle="--", alpha=0.7, zorder=1)
    plt.grid(axis="x", linestyle="--", alpha=0.7, zorder=1)
    #plt.grid(axis="x", visible=False)
    plt.gca().spines['left'].set_visible(True)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['bottom'].set_linewidth(1)
    plt.gca().spines['left'].set_linewidth(1)
    plt.gca().spines['bottom'].set_color("black")
    plt.gca().spines['left'].set_color("black")
    plt.xlim([0, 8])
    plt.ylim([0, 8])
    plt.xlabel("True cellular copy number", fontsize=20)
    plt.ylabel("Amira cellular copy number estimate", fontsize=20)
    # Adjust layout and save the file
    plt.tight_layout()
    plt.savefig(output_file, dpi=600, format='png')
    plt.savefig(output_file.replace(".png", ".pdf"), format='pdf')

def old_plot_allele_call_recall(amira_allele_results, truth_allele_names, output_file):
    rows = []

    for tech in ["r9", "r10"]:
        truth_tech = truth_allele_names.get(tech, {})
        amira_tech = amira_allele_results.get(tech, {})

        # Only compare matched samples
        matched_samples = sorted(set(truth_tech.keys()) & set(amira_tech.keys()))

        for sample in matched_samples:
            truth_counts = truth_tech.get(sample, {}) or {}
            amira_counts = amira_tech.get(sample, {}) or {}

            truth_set = {a for a, c in truth_counts.items() if c and c > 0}
            amira_set = {a for a, c in amira_counts.items() if c and c > 0}

            if len(truth_set) == 0:
                # If truth has no alleles, recall isn't meaningful; skip
                continue

            recall = len(truth_set & amira_set) / len(truth_set)

            rows.append({
                "Technology": "R9" if tech == "r9" else "R10",
                "Sample": sample,
                "Recall": float(recall),
                "Truth Alleles": len(truth_set),
                "Called Alleles": len(amira_set),
                "TP Alleles": len(truth_set & amira_set),
            })

    df = pd.DataFrame(rows)
    if df.empty:
        print("No matched samples with non-empty truth allele sets; nothing to plot.")
        return

    # Print summary (mean across samples)
    for tech_label in ["R9", "R10"]:
        vals = df.loc[df["Technology"] == tech_label, "Recall"].values
        if len(vals) > 0:
            print(f"{tech_label} mean allele-name recall (per-sample) = {np.mean(vals):.4f} (n={len(vals)})")
        else:
            print(f"{tech_label} mean allele-name recall (per-sample) = NA (n=0)")

    # ---- Plot (match the style used in your violin plot) ----
    plt.figure(figsize=(12, 12), dpi=600)
    ax = plt.gca()

    # Keep a consistent order
    tech_order = ["R9", "R10"]
    present_tech = [t for t in tech_order if t in set(df["Technology"])]

    sns.violinplot(
        y="Technology",
        x="Recall",
        data=df,
        order=present_tech,
        inner="quartile",
        linewidth=1.2,
        palette="colorblind",
        cut=0,
        bw_adjust=0.9,
        ax=ax
    )

    sns.stripplot(
        y="Technology",
        x="Recall",
        data=df,
        order=present_tech,
        jitter=0.18,
        marker="o",
        edgecolor="black",
        linewidth=0.6,
        alpha=0.65,
        size=5.5,
        palette="colorblind",
        ax=ax
    )

    # Axis styling (0 to 1)
    ax.set_xlim(0.0, 1.0)

    # Reference lines similar to your other plot
    ax.axvline(1.0, linestyle="--", linewidth=1.0, alpha=0.6)
    for v in [0.90, 0.95, 0.98, 0.99]:
        ax.axvline(v, linestyle=":", linewidth=0.8, alpha=0.30)

    ticks = sorted(set([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 0.90, 0.95, 0.98, 0.99]))
    ax.set_xticks(ticks)

    ax.set_xlabel("Allele-name Recall", fontsize=18, labelpad=10)
    ax.set_ylabel("", fontsize=18)

    ax.tick_params(axis="x", labelsize=14)
    ax.tick_params(axis="y", labelsize=14)

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_linewidth(1.0)
    ax.spines["bottom"].set_linewidth(1.0)

    ax.grid(axis="x", linestyle="--", alpha=0.30)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter, defaultdict

def plot_allele_call_recall(amira_allele_results, truth_allele_names, output_file):
    """
    Allele-name recall for Amira vs truth, split by r9 and r10.

    - apply_rules(allele_name) -> gene_type (broad type)
    - Allele-name recall is computed within gene types, but only for gene types present in BOTH truth and Amira.
    - Matching is WITHOUT replacement (Counter intersection) so duplicates are handled correctly.

    Expected input per sample:
      amira_allele_results[tech][sample] = {allele_name: count, ...}
      truth_allele_names[tech][sample]  = {allele_name: count, ...}
    """

    def group_by_gene_type(sample_counts):
        """
        sample_counts: {allele_name: count}
        returns: dict gene_type -> Counter(allele_name -> count)
        """
        gene_to_alleles = defaultdict(Counter)
        if not sample_counts:
            return gene_to_alleles

        for allele_name, cnt in sample_counts.items():
            if not (cnt and cnt > 0):
                continue
            gene_type = apply_rules(allele_name)  # broad type bucket
            gene_to_alleles[gene_type][allele_name] += int(cnt)

        return gene_to_alleles

    rows = []

    for tech in ["r9", "r10"]:
        truth_tech = truth_allele_names.get(tech, {}) or {}
        amira_tech = amira_allele_results.get(tech, {}) or {}

        matched_samples = sorted(set(truth_tech.keys()) & set(amira_tech.keys()))

        for sample in matched_samples:
            truth_counts = truth_tech.get(sample, {}) or {}
            amira_counts = amira_tech.get(sample, {}) or {}

            truth_by_gene = group_by_gene_type(truth_counts)
            amira_by_gene = group_by_gene_type(amira_counts)

            shared_gene_types = sorted(set(truth_by_gene.keys()) & set(amira_by_gene.keys()))
            if not shared_gene_types:
                continue  # you requested recall only on genes present in both

            matched_total = 0
            truth_total = 0

            for gene_type in shared_gene_types:
                truth_ctr = truth_by_gene[gene_type]
                amira_ctr = amira_by_gene[gene_type]

                # Match allele NAMES without replacement
                matched_ctr = truth_ctr & amira_ctr
                matched_total += sum(matched_ctr.values())
                truth_total += sum(truth_ctr.values())

            if truth_total == 0:
                continue

            recall = matched_total / truth_total

            rows.append({
                "Technology": "R9.4.1" if tech == "r9" else "R10.4.1",
                "Sample": sample,
                "Recall": 100*float(recall),
                "Shared_Gene_Types": len(shared_gene_types),
                "Truth_Allele_Calls_SharedGenes": int(truth_total),
                "Matched_Allele_Calls_SharedGenes": int(matched_total),
            })

    df = pd.DataFrame(rows)
    if df.empty:
        print("No matched samples with shared gene types and non-empty truth allele calls; nothing to plot.")
        return

    # Debug summary so you can sanity check it isn't artificially stuck at 1.0
    print("\nAllele-name recall (shared gene types only):")
    for tech_label in ["R9.4.1", "R10.4.1"]:
        sub = df[df["Technology"] == tech_label]
        if len(sub) == 0:
            print(f"  {tech_label}: n=0")
            continue
        print(
            f"  {tech_label}: mean={sub['Recall'].mean():.4f}, "
            f"min={sub['Recall'].min():.4f}, max={sub['Recall'].max():.4f}, n={len(sub)}"
        )

    # ---- Plot: mean bars + individual points (consistent with your style) ----
    plt.figure(figsize=(12, 12), dpi=600)
    ax = plt.gca()

    tech_order = ["R9.4.1", "R10.4.1"]
    present_tech = [t for t in tech_order if t in set(df["Technology"])]

#    sns.barplot(
#        x="Technology",
#        y="Recall",
#        data=df,
#        order=present_tech,
#        estimator=np.mean,
#        errorbar=("ci", 95),
#        palette="colorblind",
#        ax=ax,
#        zorder=2
#    )

    sns.stripplot(
        x="Technology",
        y="Recall",
        data=df,
        order=present_tech,
        jitter=0.15,
        size=12,
        palette="colorblind",
        edgecolor="black",
        linewidth=0.6,
        alpha=1,
        ax=ax,
        zorder=2
    )

    ax.set_ylim(0.0, 102)
    ax.set_ylabel("Called genes with correct allele name (%)", fontsize=20)
    ax.set_xlabel("", fontsize=20)
    ax.tick_params(axis="x", labelsize=20)
    ax.tick_params(axis="y", labelsize=20)

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_linewidth(1.0)
    ax.spines["bottom"].set_linewidth(1.0)

    ax.grid(axis="y", linestyle="--", alpha=0.7, zorder=1)

    plt.tight_layout()
    plt.savefig(output_file, dpi=900)
    if output_file.lower().endswith(".png"):
        plt.savefig(output_file[:-4] + ".pdf")
    else:
        plt.savefig(output_file + ".pdf")
    plt.close()

def old2_plot_allele_call_recall(amira_allele_results, truth_allele_names, output_file):
    """
    Per-sample allele-name recall for Amira vs truth, split by r9 and r10,
    but calculated ONLY on genes present in BOTH Amira and truth.

    Also supports duplicates (gene type and allele name can occur >1x) by doing
    matching WITHOUT replacement via Counter intersections.

    Assumptions (based on your current structures):
      - truth_allele_names[tech][sample]  is a dict-like structure where keys are gene identifiers
        (or records that can be apply_rules()'d into a gene name) and values are counts (>0 means present).
      - amira_allele_results[tech][sample] is analogous.
      - apply_rules(x) normalizes a gene/allele label for matching.
      - We treat "present" items as a multiset of labels expanded by their counts, and match with Counter.
    """

    rows = []

    for tech in ["r9", "r10"]:
        truth_tech = truth_allele_names.get(tech, {}) or {}
        amira_tech = amira_allele_results.get(tech, {}) or {}

        matched_samples = sorted(set(truth_tech.keys()) & set(amira_tech.keys()))

        for sample in matched_samples:
            truth_counts = truth_tech.get(sample, {}) or {}
            amira_counts = amira_tech.get(sample, {}) or {}

            # --- Build multiset of gene labels (expanded by counts) ---
            # (If your dict values are always 1, this still works.)
            truth_gene_multiset = []
            for k, c in truth_counts.items():
                if c and c > 0:
                    truth_gene_multiset.extend([apply_rules(k)] * int(c))

            amira_gene_multiset = []
            for k, c in amira_counts.items():
                if c and c > 0:
                    amira_gene_multiset.extend([apply_rules(k)] * int(c))

            # If either is empty, can't compute recall meaningfully
            if len(truth_gene_multiset) == 0 or len(amira_gene_multiset) == 0:
                continue

            truth_gene_counter = Counter(truth_gene_multiset)
            amira_gene_counter = Counter(amira_gene_multiset)

            # Genes present in BOTH, matched without replacement
            common_gene_counter = truth_gene_counter & amira_gene_counter  # min counts per gene label

            # If no shared genes, skip (you requested recall only on genes present in both)
            n_common_genes = sum(common_gene_counter.values())
            if n_common_genes == 0:
                continue

            # --- Build multiset of allele labels but restricted to shared genes ---
            # Your input appears to key by "allele name" already, but you said gene type + allele name can repeat.
            # Given current structure, we can only match on the normalized label (apply_rules(key)).
            # We restrict to alleles whose normalized label is in the shared-gene set, and match counts without replacement.
            shared_gene_labels = set(common_gene_counter.keys())

            truth_allele_multiset = []
            for k, c in truth_counts.items():
                if c and c > 0:
                    lab = apply_rules(k)
                    if lab in shared_gene_labels:
                        truth_allele_multiset.extend([lab] * int(c))

            amira_allele_multiset = []
            for k, c in amira_counts.items():
                if c and c > 0:
                    lab = apply_rules(k)
                    if lab in shared_gene_labels:
                        amira_allele_multiset.extend([lab] * int(c))

            truth_allele_counter = Counter(truth_allele_multiset)
            amira_allele_counter = Counter(amira_allele_multiset)

            # Match allele names without replacement
            matched_alleles_counter = truth_allele_counter & amira_allele_counter
            matched_alleles = sum(matched_alleles_counter.values())
            truth_total_restricted = sum(truth_allele_counter.values())

            if truth_total_restricted == 0:
                continue

            recall = matched_alleles / truth_total_restricted

            rows.append({
                "Technology": "R9" if tech == "r9" else "R10",
                "Sample": sample,
                "Recall": float(recall),
                "Truth_Total_Restricted": int(truth_total_restricted),
                "Matched_Alleles": int(matched_alleles),
                "Common_Genes": int(n_common_genes),
            })

    df = pd.DataFrame(rows)
    if df.empty:
        print("No matched samples with shared genes and non-empty restricted truth allele sets; nothing to plot.")
        return

    # Print summary means
    for tech_label in ["R9", "R10"]:
        vals = df.loc[df["Technology"] == tech_label, "Recall"].values
        if len(vals) > 0:
            print(f"{tech_label} mean allele-name recall (restricted to shared genes) = {np.mean(vals):.4f} (n={len(vals)})")
        else:
            print(f"{tech_label} mean allele-name recall (restricted to shared genes) = NA (n=0)")

    # ---- Plot: mean bars + individual points (match your style) ----
    plt.figure(figsize=(12, 12), dpi=600)
    ax = plt.gca()

    tech_order = ["R9", "R10"]
    present_tech = [t for t in tech_order if t in set(df["Technology"])]

    sns.barplot(
        y="Technology",
        x="Recall",
        data=df,
        order=present_tech,
        estimator=np.mean,
        errorbar=None,
        palette="colorblind",
        ax=ax
    )

    sns.stripplot(
        y="Technology",
        x="Recall",
        data=df,
        order=present_tech,
        hue="Technology",
        palette="colorblind",
        jitter=0.15,
        size=10,
        edgecolor="black",
        linewidth=0.6,
        dodge=False,
        alpha=0.8,
        ax=ax
    )

    # Remove legend entirely
    if ax.get_legend() is not None:
        ax.get_legend().remove()

    # Axis styling
    ax.set_xlim(0.0, 1.0)
    ax.set_xlabel("Allele-name Recall (shared genes only)", fontsize=20, labelpad=15)
    ax.set_ylabel("", fontsize=20)
    ax.tick_params(axis="x", labelsize=20)
    ax.tick_params(axis="y", labelsize=20)

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_linewidth(1.0)
    ax.spines["bottom"].set_linewidth(1.0)

    ax.grid(axis="x", linestyle=":", alpha=0.25)

    plt.tight_layout()
    plt.savefig(output_file, dpi=900)
    if output_file.lower().endswith(".png"):
        plt.savefig(output_file[:-4] + ".pdf")
    else:
        plt.savefig(output_file + ".pdf")
    plt.close()

def Old1_plot_allele_call_recall(amira_allele_results, truth_allele_names, output_file):
    rows = []

    for tech in ["r9", "r10"]:
        truth_tech = truth_allele_names.get(tech, {})
        amira_tech = amira_allele_results.get(tech, {})

        matched_samples = sorted(set(truth_tech.keys()) & set(amira_tech.keys()))

        for sample in matched_samples:
            truth_counts = truth_tech.get(sample, {}) or {}
            amira_counts = amira_tech.get(sample, {}) or {}
            true_gene_names = [apply_rules(g) for g in list(truth_counts.keys())]
            amira_gene_names = [apply_rules(g) for g in list(amira_counts.keys())]
            truth_set = {a for a, c in truth_counts.items() if c and c > 0}
            amira_set = {a for a, c in amira_counts.items() if c and c > 0}

            if len(truth_set) == 0:
                continue

            recall = len(truth_set & amira_set) / len(truth_set)

            rows.append({
                "Technology": "R9" if tech == "r9" else "R10",
                "Sample": sample,
                "Recall": float(recall)
            })

    df = pd.DataFrame(rows)
    if df.empty:
        print("No matched samples with non-empty truth allele sets; nothing to plot.")
        return

    # Print summary means
    for tech_label in ["R9", "R10"]:
        vals = df.loc[df["Technology"] == tech_label, "Recall"].values
        if len(vals) > 0:
            print(f"{tech_label} mean allele-name recall (per-sample) = {np.mean(vals):.4f} (n={len(vals)})")
        else:
            print(f"{tech_label} mean allele-name recall (per-sample) = NA (n=0)")

    # ---- Plot: mean bars + individual points ----
    plt.figure(figsize=(12, 12), dpi=600)
    ax = plt.gca()

    tech_order = ["R9", "R10"]
    present_tech = [t for t in tech_order if t in set(df["Technology"])]

    sns.barplot(
        y="Technology",
        x="Recall",
        data=df,
        order=present_tech,
        estimator=np.mean,
        errorbar=None,
        palette="colorblind",
        ax=ax
    )

    sns.stripplot(
        y="Technology",
        x="Recall",
        data=df,
        order=present_tech,
        hue="Technology",
        palette="colorblind",
        jitter=0.15,
        size=10,
        edgecolor="black",
        linewidth=0.6,
        dodge=False,
        ax=ax
    )

    # Remove legend entirely
    if ax.get_legend() is not None:
        ax.get_legend().remove()

    # Axis styling
    ax.set_xlim(0.0, 1.0)
    ax.set_xlabel("Allele-name Recall", fontsize=20, labelpad=15)
    ax.set_ylabel("", fontsize=20)
    ax.tick_params(axis="x", labelsize=20)
    ax.tick_params(axis="y", labelsize=20)

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_linewidth(1.0)
    ax.spines["bottom"].set_linewidth(1.0)

    ax.grid(axis="x", linestyle=":", alpha=0.25)

    plt.tight_layout()
    plt.savefig(output_file, dpi=900)
    if output_file.lower().endswith(".png"):
        plt.savefig(output_file[:-4] + ".pdf")
    else:
        plt.savefig(output_file + ".pdf")
    plt.close()

truth_dir = "truth_jsons"
amira_dir = "amira_output"
flye_dir = "AMR_finder_plus_results.flye_v2.9.3_nanopore_only_assemblies"
raven_dir = "AMR_finder_plus_results.raven_v1.8.3_nanopore_only_assemblies"
unicycler_dir = "AMR_finder_plus_results.unicycler_v0.5.0_hybrid_assemblies"
resfinder_dir = "resfinder_results"
shovill_dir = "AMR_finder_plus_results.shovill_v1.1.0_illumina_only_assemblies"
autocycler_dir = "AMR_finder_plus_results.autocycler_v0.5.0_nanopore_only_assemblies"
wtdbg2_dir = "AMR_finder_plus_results.wtdbg_v2.5_nanopore_only_assemblies"
allele_file = "/hps/nobackup/iqbal/dander/Amira_panRG_pipeline_test/Escherichia_coli_panRG_thesis/AMR_alleles_unified.fa" #"AMR_alleles_unified.fa"
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
truth_results, truth_allele_names = merge_results([f for f in glob.glob(os.path.join(truth_dir, "*.json")) if ".allele_names.json" not in f])
flye_results = process_AMRFP_results(glob.glob(os.path.join(flye_dir, "*", "*.tsv")), reference_genes)
raven_results = process_AMRFP_results(glob.glob(os.path.join(raven_dir, "*", "*.tsv")), reference_genes)
unicycler_results = process_AMRFP_results(glob.glob(os.path.join(unicycler_dir, "*", "*.tsv")), reference_genes)
resfinder_results = process_resfinder_results(glob.glob(os.path.join(resfinder_dir, "*", "ResFinder_results_tab.txt")), reference_genes)
shovill_results = process_AMRFP_results(glob.glob(os.path.join(shovill_dir, "*", "*.tsv")), reference_genes)
wtdbg_results = process_AMRFP_results(glob.glob(os.path.join(wtdbg2_dir, "*", "*.tsv")), reference_genes)
autocycler_results = process_AMRFP_results(glob.glob(os.path.join(autocycler_dir, "*", "*.tsv")), reference_genes)
# process the amira results
amira_results = {"r9": {}, "r10": {}}
amira_allele_results = {"r9": {}, "r10": {}}
samples = []
for s in glob.glob(os.path.join(amira_dir, "*")):
    if os.path.exists(os.path.join(s, "amira_results.tsv")):
        sample = os.path.basename(s)
        samples.append(sample)
        amira_table = pd.read_csv(os.path.join(s, "amira_results.tsv"), sep="\t")
        if "AUSM" in sample:
            amira_results["r10"][sample] = {}
            amira_allele_results["r10"][sample] = {}
            r10 = True
        else:
            amira_results["r9"][sample] = {}
            amira_allele_results["r9"][sample] = {}
            r10 = False
        for index, row in amira_table.iterrows():
            reference_gene = apply_rules(row["Determinant name"])
            if r10 is True:
                if reference_gene not in amira_results["r10"][sample]:
                    amira_results["r10"][sample][reference_gene] = 0
                if row["Determinant name"] not in amira_allele_results["r10"][sample]:
                    amira_allele_results["r10"][sample][row["Determinant name"]] = 0
                amira_results["r10"][sample][reference_gene] += 1
                amira_allele_results["r10"][sample][row["Determinant name"]] += 1
            if r10 is False:
                if reference_gene not in amira_results["r9"][sample]:
                    amira_results["r9"][sample][reference_gene] = 0
                if row["Determinant name"] not in amira_allele_results["r9"][sample]:
                    amira_allele_results["r9"][sample][row["Determinant name"]] = 0
                amira_results["r9"][sample][reference_gene] += 1
                amira_allele_results["r9"][sample][row["Determinant name"]] += 1
# compensate for structural variants that we are going to ignore
if "AUSMDU00021208" in amira_results["r10"]:
    if apply_rules("blaTEM-1") in amira_results["r10"]["AUSMDU00021208"]:
        amira_results["r10"]["AUSMDU00021208"][apply_rules("blaTEM-1")] = amira_results["r10"]["AUSMDU00021208"][apply_rules("blaTEM-1")] - 1

plot_allele_call_recall(amira_allele_results, truth_allele_names, os.path.join(output_dir, "figure_4e.png"))
# plot the recall and precisions of each tool
plot_recall_and_precision(truth_results,
                    [
                        amira_results,
                        flye_results,
                        raven_results,
			wtdbg_results,
                        resfinder_results,
 #                       unicycler_results,
        		autocycler_results
	          ],
                    os.path.join(output_dir, "figure_4a.png"))
# plot the recall of each single and multi copy AMR gene
plot_cn_heatmap(truth_results,
                    [
                        amira_results,
#                        unicycler_results,
                        flye_results,
                        raven_results,
			wtdbg_results,
                        resfinder_results,
			autocycler_results
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
        method = "10.4.1"
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
        amira_copy_numbers[gene_name].append(row['Approximate cellular copy number'])
    # get the nucleotide accuracy of the amira alleles
    for gene_name in true_nucleotide_sequences:
        if gene_name in amira_nucleotide_sequences:
            all_sequences = "\n".join(true_nucleotide_sequences[gene_name] + amira_nucleotide_sequences[gene_name])
            #amira_similarities = calculate_allele_accuracy_with_mafft(all_sequences, output_dir)
            amira_similarities, copy_number_tuples = calculate_allele_accuracy_with_mafft(all_sequences, output_dir, true_copy_numbers[gene_name], amira_copy_numbers[gene_name])
            all_similarities[method] += amira_similarities
            all_copy_number_tuples[method] += copy_number_tuples

#plot_nucleotide_results(all_similarities, os.path.join(output_dir, "nucleotide_accuracies.png"))
plot_nucleotide_results_violin({
                            "Flye AMRFP": calculate_assembler_accuracies("flye_v2.9.3_nanopore_only"),
                            "Raven AMRFP": calculate_assembler_accuracies("raven_v1.8.3_nanopore_only"),
                            "Wtdbg2 AMRFP": calculate_assembler_accuracies("wtdbg_v2.5_nanopore_only"),
                            "Autocycler AMRFP": calculate_assembler_accuracies("autocycler_v0.5.0_nanopore_only"),
			    "Amira": all_similarities},
                            os.path.join(output_dir, "figure_4c.png"))
plot_copy_numbers(all_copy_number_tuples, os.path.join(output_dir, "figure_4b.png"))

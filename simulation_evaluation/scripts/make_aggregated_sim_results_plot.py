import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import statistics
import subprocess
from Bio import AlignIO
from tqdm import tqdm
import numpy as np
import json

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
    if "oqx" in gene:
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

def process_AMRFP_results(file_list):
    merged_results = {}
    for f in file_list:
        sample = os.path.basename(os.path.dirname(f))
        merged_results[sample] = {}
        results = pd.read_csv(f, sep="\t")
        for index, row in results.iterrows():
            gene = row["Gene symbol"]
            if row["Element subtype"] == "POINT":
                continue
            gene = apply_rules(gene)
            if gene not in merged_results[sample]:
                merged_results[sample][gene] = 0
            merged_results[sample][gene] += 1
    return merged_results[sample]

def calculate_gene_accuracy(truth, amira):
    tp = min(truth, amira)
    fn = max(0, truth - amira)
    fp = max(0, amira-truth)
    return tp / (tp + fn)# + fp)

def calculate_gene_precision(truth, amira):
    tp = min(truth, amira)
    fp = max(0, amira-truth)
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    if amira != 0:
        return precision
    else:
        return Nones

def plot_sim_recall_results(data_list, titles, x_label, y_label, output_plot):
    # Create a 2x3 grid for the plots
    plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 30})
    fig, axes = plt.subplots(6, 1, figsize=(22, 45))
    axes = axes.flatten()

    # Line styles and marker styles for different read lengths
    line_styles = ['-', '--', '-.', ':']  # Line styles for the lines
    marker_styles = ['o', 's', '^', 'D']  # Markers for the points

    # Create an empty list to hold lines and labels for the global legend
    lines = []
    labels = []

    # Loop through each scenario and plot the recall results for each read length
    for i, s in enumerate(data_list.keys()):
        for j, l in enumerate(lengths):
            depth_data = [data_list[s][d][j] for d in depths]  # Recall data for this read length
            # Plot with both line style and marker
            line, = axes[i].plot(depths, depth_data, line_styles[j], marker=marker_styles[j], markersize=14, label=f"{l}kb")
            if i == 0:  # Only collect lines and labels from the first subplot
                lines.append(line)
                labels.append(f"{l}kb")
        axes[i].set_title(titles[i], fontsize=30)
        axes[i].set_ylim([0, 1.05])

        # Remove the top and right spines (box around the plot)
        axes[i].spines['top'].set_visible(False)
        axes[i].spines['right'].set_visible(False)

        # Increase the thickness of the bottom and left axis lines
        axes[i].spines['bottom'].set_linewidth(2)
        axes[i].spines['left'].set_linewidth(2)

        # Increase the thickness of the ticks
        axes[i].tick_params(axis='both', width=2)

        # Add y-axis grid with specified alpha and linewidth
        axes[i].grid(axis='y', linestyle='-', alpha=0.3, linewidth=1)
        if i == len(data_list) - 1:
            axes[i].set_xlabel(x_label, fontsize=30)

    # Add unified x and y labels for the entire figure, bringing them closer
    #fig.supxlabel(x_label, fontsize=30, x=0.5, y=0.06)  # Adjust y to move the x label closer
    fig.supylabel(y_label, fontsize=30, x=0.04, y=0.5)  # Adjust x to move the y label closer

    # Adjust the layout to make room for the labels and legend
    plt.subplots_adjust(left=0.12, right=0.85, top=0.95, bottom=0.12, wspace=0.3, hspace=0.5)

    # Create a single legend outside the subplots to the right without the box (frame)
    legend = fig.legend(lines, labels, loc='center right', ncol=1, title="Mean read length", frameon=False, bbox_to_anchor=(1.05, 0.5))

    # Save the figure with the legend and labels included
    plt.savefig(output_plot, dpi=600, bbox_inches='tight')
    plt.savefig(output_plot.replace(".png", ".pdf"), bbox_inches='tight')


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

def plot_sim_nucleotide_results(similarity_dict, output_file):
    # Convert the dictionary into a DataFrame for easier plotting
    data = []
    for depth, similarities in similarity_dict.items():
        print(depth, statistics.mean(similarities))
        for similarity in similarities:
            data.append({'Depth': depth, 'Similarity': similarity})
    df = pd.DataFrame(data)

    # Set font style and size
    plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 16})

    # Create the plot
    fig, ax = plt.subplots(figsize=(12, 6))

    # Create the box plot
    sns.boxplot(
        x='Depth',
        y='Similarity',
        data=df,
        linewidth=1.5,
        palette="colorblind"  # Using a muted color palette for a cleaner look
    )

    # Set title and axis labels
    ax.set_ylabel('Allele accuracy', labelpad=10)
    ax.set_xlabel('Read depth', labelpad=10)

    # Set y-axis limits
    #ax.set_ylim([0.95, 1])

    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Increase the thickness of the bottom and left axis lines
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)

    # Increase the thickness of the ticks
    ax.tick_params(axis='both', width=2)

    # Add a y-axis grid with specified alpha and linewidth
    ax.grid(axis='y', linestyle='-', alpha=0.3, linewidth=1)

    # Clean up the grid lines and adjust layout
    plt.tight_layout()

    # Save the plot to the output file
    plt.savefig(output_file, dpi=600, bbox_inches='tight')
    plt.savefig(output_file.replace(".png", ".pdf"), bbox_inches='tight')

def plot_combined_results(data_list, titles, x_label, y_label, nucleotide_similarity_dict, output_plot):
    # Create a 3x2 grid for the recall plots and 1x1 for nucleotide accuracy (combined into one grid)
    plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 16})
    fig = plt.figure(figsize=(12, 20))  # Adjusted figure size for additional plot

    # Create a grid for subplots, 4 rows and 2 columns (3x2 for recall plots, 1x2 for nucleotide accuracy)
    grid = fig.add_gridspec(4, 2, height_ratios=[1, 1, 1, 0.5])  # Last row smaller for nucleotide accuracy

    axes = [fig.add_subplot(grid[i, j]) for i in range(3) for j in range(2)]  # Axes for recall plots
    nucleotide_ax = fig.add_subplot(grid[3, :])  # Axes for nucleotide accuracy

    # Line styles and marker styles for different read lengths
    line_styles = ['-', '--', '-.', ':']  # Line styles for the lines
    marker_styles = ['o', 's', '^', 'D']  # Markers for the points

    # Create an empty list to hold lines and labels for the global legend
    lines = []
    labels = []

    # Loop through each scenario and plot the recall results for each read length
    for i, s in enumerate(data_list.keys()):
        for j, l in enumerate(lengths):
            depth_data = [data_list[s][d][j] for d in depths]  # Recall data for this read length
            # Plot with both line style and marker
            line, = axes[i].plot(depths, depth_data, line_styles[j], marker=marker_styles[j], label=f"{l}kb")
            if i == 0:  # Only collect lines and labels from the first subplot
                lines.append(line)
                labels.append(f"{l}kb")
        axes[i].set_title(titles[i], fontsize=16)
        axes[i].set_ylim([0, 1.05])

        # Remove the top and right spines (box around the plot)
        axes[i].spines['top'].set_visible(False)
        axes[i].spines['right'].set_visible(False)

        # Increase the thickness of the bottom and left axis lines
        axes[i].spines['bottom'].set_linewidth(2)
        axes[i].spines['left'].set_linewidth(2)

        # Increase the thickness of the ticks
        axes[i].tick_params(axis='both', width=2)

        # Add y-axis grid with specified alpha and linewidth
        axes[i].grid(axis='y', linestyle='-', alpha=0.3, linewidth=1)

    # Plot the nucleotide accuracy in the bottom section
    # Convert the dictionary into a DataFrame for easier plotting
    data = []
    for depth, similarities in nucleotide_similarity_dict.items():
        for similarity in similarities:
            data.append({'Depth': depth, 'Similarity': similarity})
    df = pd.DataFrame(data)

    # Create the box plot for nucleotide accuracy
    sns.boxplot(
        x='Depth',
        y='Similarity',
        data=df,
        linewidth=1.5,
        palette="colorblind",
        ax=nucleotide_ax  # Plot on the nucleotide accuracy subplot
    )

    # Set axis labels and title for nucleotide accuracy
    nucleotide_ax.set_ylabel('Allele accuracy', labelpad=10)
    nucleotide_ax.set_xlabel('Read depth', labelpad=10)
    nucleotide_ax.set_ylim([0.95, 1])
    nucleotide_ax.spines['top'].set_visible(False)
    nucleotide_ax.spines['right'].set_visible(False)
    nucleotide_ax.spines['bottom'].set_linewidth(2)
    nucleotide_ax.spines['left'].set_linewidth(2)
    nucleotide_ax.tick_params(axis='both', width=2)
    nucleotide_ax.grid(axis='y', linestyle='-', alpha=0.3, linewidth=1)

    # Add unified x and y labels for the entire figure, bringing them closer
    fig.supxlabel(x_label, fontsize=16, x=0.5, y=0.04)
    #fig.supylabel(y_label, fontsize=16, x=0.04, y=0.5)

    # Adjust the layout to make room for the labels and legend
    plt.subplots_adjust(left=0.12, right=0.85, top=0.95, bottom=0.12, wspace=0.3, hspace=0.3)

    # Create a single legend outside the subplots to the right without the box (frame)
    legend = fig.legend(lines, labels, loc='center right', ncol=1, title="Mean read length", bbox_to_anchor=(1, 0.5), frameon=False)

    # Save the figure with the legend and labels included
    plt.savefig(output_plot, dpi=600, bbox_inches='tight')
    plt.savefig(output_plot.replace(".png", ".pdf"), bbox_inches='tight')

def plot_copy_numbers(copy_number_tuples_by_depth, output_file):
    plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 16})
    fig = plt.figure(figsize=(12, 12))  # Adjusted figure size for additional plot
    line_styles = ['-', '--', '-.', ':']  # Line styles for the lines
    marker_styles = ['o', 's', '^', 'D']  # Markers for the points
    for i, d in enumerate(copy_number_tuples_by_depth):
        x_vals = [float(t[0]) for t in copy_number_tuples_by_depth[d]]
        y_vals = [float(t[1]) for t in copy_number_tuples_by_depth[d]]
        # Scatter plot for data points
        plt.scatter(x_vals, y_vals, label=f"{d}x", marker=marker_styles[i % len(marker_styles)])
        # Calculate line of best fit
        if len(x_vals) > 1:  # Ensure there's more than one point to fit a line
            m, b = np.polyfit(x_vals, y_vals, 1)  # Linear regression (slope, intercept)
            best_fit_line = [m * x + b for x in x_vals]
            # Plot the line of best fit
            plt.plot(x_vals, best_fit_line, linestyle=line_styles[i % len(line_styles)], label=f"{d}x fit")
    plt.xlim([0, 6])
    plt.ylim([0, 6])
    plt.xlabel("True contig copy number")
    plt.ylabel("Amira contig copy number")  # Corrected to ylabel
    plt.legend()  # Added legend to distinguish depths and their best-fit lines
    plt.savefig(output_file, dpi=600, bbox_inches='tight')
    plt.savefig(output_file.replace(".png", ".pdf"), bbox_inches='tight')

def process_AMRFP_results(f, reference_genes):
    AMRFP_results = {}
    sample = os.path.basename(os.path.dirname(f))
    results = pd.read_csv(f, sep="\t")
    for index, row in results.iterrows():
        gene = row["Gene symbol"]
        if row["Element subtype"] == "POINT":
            continue
        gene = apply_rules(gene)
        if gene not in reference_genes:
            continue
        if gene not in AMRFP_results:
            AMRFP_results[gene] = 0
        AMRFP_results[gene] += 1
    return AMRFP_results

output_dir = "aggregated_sim_results"
# load the reference AMR genes
with open("AMR_alleles_unified.fa") as i:
    allele_rows = i.read().split(">")[1:]
reference_genes = set()
for r in allele_rows:
    if r != "":
        amira_allele, reference_allele = r.split("\n")[0].split(";")
        reference_genes.add(apply_rules(reference_allele.split(".NG")[0]))
# make the output dir
if not os.path.exists(output_dir):
    os.mkdir(output_dir)
# iterate through the different depths and lengths for this sample
depths = [5, 10, 20, 40, 80]
lengths = [5, 10, 20, 40]
scenario = [0, 1, 2, 4, 3, 5]

# Initialize dictionary to store data per scenario, depth, and length
data_list = {s: {d: [] for d in depths} for s in scenario}
data_list_flye = {s: {d: [] for d in depths} for s in scenario}

def rc_sequence(sequence):
    replacement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(list(reversed([replacement[b] for b in list(sequence)])))

all_similarities = []
all_similarities_by_depth = {d:[] for d in depths}
copy_number_tuples_by_depth = {d:[] for d in depths}
for s in tqdm(scenario):
    # Import the truth tsv
    truth_df = pd.read_csv(os.path.join("simulated_assemblies", f"test_{s}.AMR_genes.tsv"), sep="\t")
    truth_counts = {}
    true_gene_positions = {}
    for index, row in truth_df.iterrows():
        gene_name = apply_rules(row["Allele name"].split(";")[1].split(".")[0])
        if gene_name not in truth_counts:
            truth_counts[gene_name] = 0
        truth_counts[gene_name] += 1
        if row["Contig"] not in true_gene_positions:
            true_gene_positions[row["Contig"]] = {}
        true_gene_positions[row["Contig"]][(row["Start"], row["End"], row["Copy number"], row["Strand"])] = gene_name
    # import the truth fasta
    with open(os.path.join("simulated_assemblies", f"test_{s}.fasta")) as i:
        chromosome, plasmid = i.read().split(">")[1:3]
    true_nucleotide_sequences = {}
    true_copy_numbers = {}
    for contig in true_gene_positions:
        for start, end, copy_number, strand in true_gene_positions[contig]:
            if contig == "chromosome":
                seq = "".join(chromosome.split("\n")[1:])
            elif contig == "plasmid":
                seq = "".join(plasmid.split("\n")[1:])
            gene_name = true_gene_positions[contig][(start, end, copy_number, strand)]
            if gene_name not in true_nucleotide_sequences:
                true_nucleotide_sequences[gene_name] = []
                true_copy_numbers[gene_name] = []
            if strand == "+":
                truth_sequence = seq[start-1: end]
            elif strand == "-":
                truth_sequence = rc_sequence(seq[start-1: end])
            true_nucleotide_sequences[gene_name].append(f">{gene_name}_truth\n{truth_sequence}")
            true_copy_numbers[gene_name].append(copy_number)
    # Get the amira results for each read length and depth
    for l in lengths:
        for d in depths:
            amira_out = os.path.join(f"simulations_{l}kb_r10", f"amira_{d}", f"test_{s}")
            flye_out = os.path.join(f"simulations_{l}kb_r10", f"AMRFP_flye_{d}", f"test_{s}")
            if not os.path.exists(os.path.join(amira_out, "amira_results.tsv")):
                data_list[s][d].append(0)
                continue
            if not os.path.exists(os.path.join(flye_out, "AMR_finder_plus_results.tsv")):
                data_list_flye[s][d].append(0)
            # load the amira data
            amira_df = pd.read_csv(os.path.join(amira_out, "amira_results.tsv"), sep="\t")
            amira_counts = {}
            amira_nucleotide_sequences = {}
            amira_copy_numbers = {}
            for index, row in amira_df.iterrows():
                gene_name = apply_rules(row["Determinant name"])
                if gene_name not in amira_counts:
                    amira_counts[gene_name] = 0
                amira_counts[gene_name] += 1
                # get the nucleotide sequences of each gene
                amira_fasta = os.path.join(amira_out, "AMR_allele_fastqs", row["Amira allele"], "06.final_sequence.fasta")
                with open(amira_fasta) as i:
                    seq = "\n".join(i.read().split("\n")[1:])
                if gene_name not in amira_nucleotide_sequences:
                    amira_nucleotide_sequences[gene_name] = []
                    amira_copy_numbers[gene_name] = []
                amira_nucleotide_sequences[gene_name].append(f">{gene_name}_amira\n{seq}")
                amira_copy_numbers[gene_name].append(row["Approximate copy number"])
            # Calculate the recall for each gene
            gene_recalls_amira = []
            for g in set(list(amira_counts.keys()) + list(truth_counts.keys())):
                if g not in truth_counts:
                    continue
                r_amira = calculate_gene_accuracy(truth_counts[g], amira_counts.get(g, 0))
                if r_amira is not None:
                    gene_recalls_amira.append(r_amira)
            try:
                mean_recall = statistics.mean(gene_recalls_amira)
            except:
                mean_recall = 1
            # Store the mean recall for this scenario, depth, and length
            data_list[s][d].append(mean_recall)
            # get the nucleotide accuracy of the amira alleles
            for gene_name in true_nucleotide_sequences:
                if gene_name in amira_nucleotide_sequences:
                    all_sequences = "\n".join(true_nucleotide_sequences[gene_name] + amira_nucleotide_sequences[gene_name])
                    amira_similarities, copy_number_tuples = calculate_allele_accuracy_with_mafft(all_sequences, output_dir, true_copy_numbers[gene_name], amira_copy_numbers[gene_name])
                    all_similarities += amira_similarities
                    all_similarities_by_depth[d] += amira_similarities
                    copy_number_tuples_by_depth[d] += copy_number_tuples
            # load the flye data
            if not os.path.exists(os.path.join(flye_out, "AMR_finder_plus_results.tsv")):
                continue
            flye_counts = process_AMRFP_results(os.path.join(flye_out, "AMR_finder_plus_results.tsv"), reference_genes)
            # Calculate the recall for each gene
            gene_recalls_flye = []
            for g in set(list(flye_counts.keys()) + list(truth_counts.keys())):
                if g not in truth_counts:
                    continue
                r_flye = calculate_gene_accuracy(truth_counts[g], flye_counts.get(g, 0))
                if r_flye is not None:
                    gene_recalls_flye.append(r_flye)
            try:
                mean_flye_recall = statistics.mean(gene_recalls_flye)
            except:
                mean_flye_recall = 1
            # Store the mean recall for this scenario, depth, and length
            data_list_flye[s][d].append(mean_flye_recall)

scenario_mapping = {0: 1, 1 : 2, 2 : 3, 4 : 4, 3 : 5, 7 : 6}
modified_data_list = {}
flye_modified_data_list = {}
for s in data_list:
    modified_data_list[scenario_mapping[s]] = data_list[s]
for s in data_list_flye:
    flye_modified_data_list[scenario_mapping[s]] = data_list_flye[s]
with open("recall_plot_inputs.json", "w") as o:
    o.write(json.dumps(modified_data_list))
with open("flye_recall_plot_inputs.json", "w") as o:
    o.write(json.dumps(flye_modified_data_list))
titles = [f"Scenario {scenario_mapping[s]}" for s in scenario]
plot_copy_numbers(copy_number_tuples_by_depth, os.path.join(output_dir, "copy_numbers_r10.png"))
plot_sim_recall_results(data_list, titles, "Read depth", "AMR gene recall", os.path.join(output_dir, "simulation_recall_r10.png"))
plot_sim_nucleotide_results(all_similarities_by_depth, os.path.join(output_dir, "simulation_nucleotide_accuracy_r10.png"))
plot_combined_results(data_list, titles, "Read depth", "AMR gene recall", all_similarities_by_depth, os.path.join(output_dir, "combined_simulation_r10.png"))
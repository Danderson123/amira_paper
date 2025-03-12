import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import pandas as pd
import numpy as np
from matplotlib.patches import FancyArrow
import json
import matplotlib.colors as mc
import matplotlib.image as image
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
import subprocess
from Bio import AlignIO
import statistics
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from palettable import cartocolors


# Define constants and shared settings
plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 30})

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

# Function to generate colorblind-safe colors
def generate_colorblind_safe_colors(n):
    colors = sns.color_palette("colorblind", n_colors=n)
    return colors

# Plot gene positions
def plot_gene_positions(ax, genes_dict, positions_dict, gene_colors, scenario, global_max_pos):
    ax.set_xlim(-500, global_max_pos)
    ax.set_ylim(0, len(genes_dict) + 1)
    read_num = 1
    fixed_head_width = 8
    fixed_tail_width = 5
    # Convert fixed arrow widths to data coordinates
    inv = ax.transData.inverted()
    fig_height_in_data_coords = inv.transform((0, fixed_head_width))[1] - inv.transform((0, 0))[1]
    for seq_id, gene_list in genes_dict.items():
        pos_list = positions_dict[seq_id]
        first_gene_start, last_gene_end = pos_list[0][0], pos_list[-1][1]
        ax.hlines(y=read_num, xmin=first_gene_start, xmax=last_gene_end, colors='black', linewidth=2, alpha=0.8)
        for gene, (start, end) in zip(gene_list, pos_list):
            gene_name = gene[1:]
            facecolor = gene_colors.get(gene_name, "gray")
            strand = gene[0]
            arrow_length = end - start
            head_length = min(600, arrow_length * 0.8)  # Adjust head length relative to arrow size
            # Adjust the arrow based on the strand (+ or -)
            if strand == '+':
                if arrow_length >= 600:
                    # Forward strand: arrow starts at 'start' and points to 'end'
                    ax.add_patch(FancyArrow(
                        start, read_num, arrow_length, 0,
                        width=fig_height_in_data_coords * fixed_tail_width, 
                        head_width=fig_height_in_data_coords * fixed_head_width, 
                        head_length=head_length,
                        length_includes_head=True, facecolor=facecolor, edgecolor='black', linewidth=2, zorder=2
                    ))
                else:
                    ax.add_patch(FancyArrow(
                        start, read_num, arrow_length, 0,
                        width=fig_height_in_data_coords * fixed_tail_width, 
                        head_width=fig_height_in_data_coords * fixed_head_width, 
                        head_length=head_length,
                        length_includes_head=True, facecolor=facecolor, edgecolor='black', linewidth=2, zorder=2
                    ))
            else:
                if arrow_length >= 600:
                    # Reverse strand: arrow starts at 'end' and points to 'start'
                    ax.add_patch(FancyArrow(
                        end, read_num, -arrow_length, 0,
                        width=fig_height_in_data_coords * fixed_tail_width, 
                        head_width=fig_height_in_data_coords * fixed_head_width, 
                        head_length=head_length,
                        length_includes_head=True, facecolor=facecolor, edgecolor='black', linewidth=2, zorder=2
                    ))
                else:
                    ax.add_patch(FancyArrow(
                        end, read_num, -arrow_length, 0,
                        width=fig_height_in_data_coords * fixed_tail_width, 
                        head_width=fig_height_in_data_coords * fixed_head_width, 
                        head_length=head_length,
                        length_includes_head=True, facecolor=facecolor, edgecolor='black', linewidth=2, zorder=2
                    ))
        read_num += 1
    # Set custom y-tick labels for the specific scenario, if available
    custom_y_labels = {
        1: [" "],
        2: ["chromosome"],
        3: ["plasmid", "chromosome"],
        4: ["chromosome"],
        5: ["plasmid", "chromosome"],
        6: ["plasmid", "chromosome", "chromosome"]
    }
    if scenario in custom_y_labels:
        ax.set_yticks(range(1, len(custom_y_labels[scenario]) + 1))
        ax.set_yticklabels(custom_y_labels[scenario])
        if scenario != 6:
            ax.set_xticklabels("")

# Plot AMR gene recall
def plot_amr_recall(ax, depths, recalls, lengths, scenario):
    line_styles = ['-', '--', '-.', ':']
    marker_styles = ['o', 's', '^', 'D']
    for j, l in enumerate(lengths):
        recall_data = [recalls[str(d)][j] for d in depths]
        ax.plot(depths, recall_data, line_styles[j], marker=marker_styles[j], markersize=14, label=f"{l}kb")
    ax.set_ylim(0, 1.05)
    if scenario != "6":
        ax.set_xticklabels("")
    ax.grid(axis='y', linestyle='-', alpha=0.3, linewidth=1)

def calculate_gene_accuracy(truth, amira):
    tp = min(truth, amira)
    fn = max(0, truth - amira)
    return tp / (tp + fn)

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

def adjust_lightness(color, amount=0.5):
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])

def plot_combined(genes_data, amira_recall_data, flye_recall_data, output_file):
    depths = [5, 10, 20, 40, 80]
    lengths = [5, 10, 20, 40]
    scenarios = list(genes_data.keys())
    fig = plt.figure(figsize=(30, 45))
    grid = gridspec.GridSpec(5, 2, hspace=0.1, wspace=0.2, width_ratios=[1, 1])
    global_max_end = 0
    for s in genes_data:
        _, positions_dict, _ = genes_data[s]
        if len(positions_dict) != 0:
            for contig in positions_dict:
                max_end = max([end for start, end in positions_dict[contig]])
                if max_end > global_max_end:
                    global_max_end = max_end
    # Function used to normalize values into 0-1 scale.
    palette = sns.color_palette("colorblind", n_colors=len(depths))
    depth_color_map = {depths[d]: palette[d] for d in range(len(depths))}
    for i, scenario in enumerate(scenarios):
        if scenario == 1:
            continue
        genes_dict, positions_dict, gene_colors = genes_data[scenario]
        recalls_amira = amira_recall_data[scenario]
        recalls_flye = flye_recall_data[scenario]
        # Context plot (left column)
        ax1 = fig.add_subplot(grid[i-1, 0])
        plot_gene_positions(ax1, genes_dict, positions_dict, gene_colors, scenario, global_max_end)
        if scenario == 6:
            ax1.set_xlabel("Nucleotide position (bp)")

        # Lollipop plot for recall differences (right column)
        ax2 = fig.add_subplot(grid[i-1, 1])
        plot_dict = {"Lengths": [], "Depths": [], "Amira": [], "Flye": [], "x": [], "color": []}
        for depth in recalls_amira:
            for i in range(len(recalls_amira[depth])):
                plot_dict["Lengths"].append(int(lengths[i]))
                plot_dict["Depths"].append(int(depth))
                plot_dict["Amira"].append(recalls_amira[depth][i])
                try:
                    plot_dict["Flye"].append(recalls_flye[depth][i])
                except:
                    plot_dict["Flye"].append(0)
                plot_dict["x"].append((depths.index(int(depth)) + 1)*25 + 200*(i+1))
                plot_dict["color"].append(depth_color_map[int(depth)])
        plot_df = pd.DataFrame(plot_dict)
        ax2.vlines(
            x="x",
            ymin="Amira",
            ymax="Flye",
            data = plot_df,
            color=list(plot_df["color"])
        )
        ax2.scatter(
            "x",
            "Amira",
            data=plot_df,
            color=list(plot_df["color"]),
            zorder=3,
            s=200,
            marker="o"
        )
        ax2.scatter(
            "x",
            "Flye",
            data=plot_df,
            color=list(plot_df["color"]),
            zorder=3,
            s=200,
            marker="x"
        )

        ax2.set_ylim([-0.05, 1.05])
        xticks_positions = [275, 475, 675, 875]
        custom_labels = [5000, 10000, 20000, 40000]
        if scenario == 6:
            ax2.set_xlabel("Mean read length (bp)")
            ax2.set_xticks(xticks_positions)
            ax2.set_xticklabels(custom_labels)
        else:
            ax2.set_xticks(xticks_positions)
            ax2.set_xticklabels([])
        if scenario == 4:
            ax2.set_ylabel("AMR gene recall")
        # Adjust aesthetics for plots
        ax1.spines['left'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['top'].set_visible(False)
        ax1.spines['bottom'].set_linewidth(True)
        ax1.grid(axis="y", visible=False)
        ax1.grid(axis="x", linestyle="--", alpha=0.7, zorder=1)
        ax2.grid(axis="y", linestyle="--", alpha=0.7, zorder=1)
        ax2.grid(axis="x", visible=False)
        ax2.spines['left'].set_visible(True)
        ax2.spines['right'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax2.spines['bottom'].set_linewidth(True)

    # Save plots
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
        if row["Method"] == "PARTIALX" or row["Method"] == "PARTIALP":
            continue
        gene = apply_rules(gene)
        if gene not in reference_genes:
            continue
        if gene not in AMRFP_results:
            AMRFP_results[gene] = 0
        AMRFP_results[gene] += 1
    return AMRFP_results

def parse_gff(file_path):
    genes_dict = {}
    positions_dict = {}
    amr_alleles = set()
    with open(file_path) as f:
        for line in f:
            if line.startswith(">"):
                break
            if line.startswith("#") or line.strip() == "":
                continue
            parts = line.strip().split('\t')
            seq_id, source, feature_type, start, end, score, strand, phase, attributes = parts

            attr_dict = {attributes.split(';')[0].split('=')[0]: attributes.split(';')[0].split('=')[1]}
            gene_name = attr_dict.get('Name', 'Unknown')
            if "AMR_alleles" in attributes:
                amr_alleles.add(gene_name)
            if seq_id not in genes_dict:
                genes_dict[seq_id] = []
                positions_dict[seq_id] = []

            # Store the current start and end positions
            start = int(start)
            end = int(end)

            # If this is the first entry for this seq_id, check if an offset is needed
            if len(positions_dict[seq_id]) == 0 and start != 0:
                offset = start  # Record the offset if the first start position isn't 0
            else:
                offset = 0  # No offset if the first start position is already 0

            # Adjust the start and end positions by subtracting the offset
            adjusted_start = start - offset
            adjusted_end = end - offset

            genes_dict[seq_id].append(strand + gene_name)
            positions_dict[seq_id].append((adjusted_start, adjusted_end))
    return genes_dict, positions_dict, amr_alleles

output_dir = "simulation_results"
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

scenario_mapping = {0: 1, 1 : 2, 2 : 3, 4 : 4, 3 : 5, 5 : 6}
modified_data_list = {}
flye_modified_data_list = {}
for s in data_list:
    modified_data_list[scenario_mapping[s]] = data_list[s]
for s in data_list_flye:
    flye_modified_data_list[scenario_mapping[s]] = data_list_flye[s]

# Directory containing GFF files
gff_directory = "context_gffs"
output_file_base = "simulation_results/combined.png"

# Load all GFF files from the directory
genes_dicts = []
positions_dicts = []
amr_alleles = set()
scenario_mapping = {"test_0": "test_0", "test_1" : "test_1", "test_2" : "test_2", "test_3" : "test_4", "test_4" : "test_3", "test_5" : "test_5"}
scenarios = []
import glob
gff_files = [f for f in glob.glob(os.path.join(gff_directory, "*.gff")) if os.path.basename(f).replace(".gff", "") in scenario_mapping]
gff_files = list(sorted(gff_files, key=lambda x: scenario_mapping[os.path.basename(x).replace(".gff", "")]))
for file_path in gff_files:
    genes_dict, positions_dict, amr = parse_gff(file_path)
    if os.path.basename(file_path) != "test_3.gff" and os.path.basename(file_path) != "test_5.gff":
        genes_dicts.append(genes_dict)
        positions_dicts.append(positions_dict)
    if os.path.basename(file_path) == "test_3.gff":
        reversed_by_keys = dict(list(reversed(list(genes_dict.items()))))
        genes_dicts.append(reversed_by_keys)
        positions_dicts.append(positions_dict)
    if os.path.basename(file_path) == "test_5.gff":
        reversed_by_keys = dict(list(reversed(list(genes_dict.items()))))
        genes_dicts.append(reversed_by_keys)
        positions_dicts.append(positions_dict)
    amr_alleles.update(amr)
    scenarios.append(scenario_mapping[os.path.basename(file_path).replace(".gff", "")])
# Generate unique colors for each gene
unique_genes = set(g[1:] for genes_dict in genes_dicts for genes in genes_dict.values() for g in genes)
gene_colors = {gene: "lightblue" if gene not in amr_alleles else "red" for gene in unique_genes}

genes_data = {}
for i in range(len(genes_dicts)):
    scenario = i + 1
    genes_data[scenario] = (genes_dicts[i], positions_dicts[i], gene_colors)

plot_combined(genes_data, modified_data_list, flye_modified_data_list, 'simulation_results/final_combined_plot.png')
for depth, similarities in all_similarities_by_depth.items():
    print(depth, statistics.mean(similarities))
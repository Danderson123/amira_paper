import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrow, FancyArrowPatch
import seaborn as sns
import os
import matplotlib.lines as mlines
import matplotlib.gridspec as gridspec
import json

def generate_colorblind_safe_colors(n):
    """Generate n colorblind-safe colors using a colorblind-friendly colormap."""
    colors = sns.color_palette("colorblind", n_colors=n)
    return colors

def plot_reads_gene_positions(genes_dicts, positions_dicts, output_file_base, scenarios, gene_colors):
    """Plot gene positions for each read in each scenario, with consistent X-axis limits across subplots,
       and the height of each subplot scales proportionally with the number of rows in the subplot, but arrow widths remain consistent."""
    
    num_scenarios = len(scenarios)
    plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 30})
    
    # Calculate the number of rows for each subplot (i.e., number of genes/reads)
    num_rows = [len(genes_dict) for genes_dict in genes_dicts]
    
    # Create a gridspec layout with heights proportional to the number of rows in each scenario
    total_height = sum([n for n in num_rows]) * 5  # Adjust the total figure height based on the number of rows
    fig = plt.figure(figsize=(22, total_height))
    gs = gridspec.GridSpec(num_scenarios, 1, height_ratios=[1 for n in num_rows], hspace=0.5)

    # Define custom y-tick labels for specific scenarios
    custom_y_labels = {
        "test_0": [],
        "test_1": ["chromosome"],
        "test_2": ["plasmid", "chromosome"],
        "test_3": ["chromosome"],
        "test_4": ["plasmid", "chromosome"],
        "test_5": ["plasmid", "chromosome", "chromosome"]
    }

    # Find the global maximum end position across all reads and scenarios
    global_max_pos = max(end for positions_dict in positions_dicts
                         for positions in positions_dict.values()
                         for _, end in positions)

    # Fixed arrow height in figure-relative coordinates
    fixed_head_width = 8 # Proportional to figure height
    fixed_tail_width = 5

    for i, (genes_dict, positions_dict, scenario) in enumerate(zip(genes_dicts, positions_dicts, scenarios)):
        ax = fig.add_subplot(gs[i])
        ax.set_title(f"Scenario: {scenario.replace('test_', '')}", fontsize=30)
        ax.set_xlim(-500, global_max_pos)  # Set the X-axis limit based on the global maximum end position
        ax.set_ylim(0, len(genes_dict) + 1)
        # Convert fixed arrow widths to data coordinates
        inv = ax.transData.inverted()
        fig_height_in_data_coords = inv.transform((0, fixed_head_width))[1] - inv.transform((0, 0))[1]
        read_num = 1
        for seq_id, gene_list in genes_dict.items():
            pos_list = positions_dict[seq_id]
            # Draw a gray line that spans the positions of the gene segments only
            first_gene_start, last_gene_end = pos_list[0][0], pos_list[-1][1]
            # Draw the gray line behind the arrows only within the span of the gene segments
            ax.hlines(y=read_num, xmin=first_gene_start, xmax=last_gene_end, colors='black', linewidth=2, alpha=0.8)
            for gene, (start, end) in zip(gene_list, pos_list):
                gene_name = gene[1:]  # Exclude the strand symbol
                facecolor = gene_colors.get(gene_name, "gray")
                strand = gene[0]
                # Calculate arrow length, subtracting a small portion to account for the arrowhead size
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
        if scenario in custom_y_labels:
            ax.set_yticks(range(1, len(custom_y_labels[scenario]) + 1))
            ax.set_yticklabels(custom_y_labels[scenario])
            if scenario != "test_5":
                ax.set_xticklabels("")

    # Add a unified x and y label using suplabels
    fig.supxlabel('Position', fontsize=30)

    plt.tight_layout()
    plt.savefig(output_file_base, dpi=600)
    plt.savefig(output_file_base.replace(".png", ".pdf"))


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

# Directory containing GFF files
gff_directory = "context_gffs"
output_file_base = "aggregated_sim_results/combined.png"

# Load all GFF files from the directory
genes_dicts = []
positions_dicts = []
amr_alleles = set()
scenario_mapping = {"test_0": "test_0", "test_1" : "test_1", "test_2" : "test_2", "test_3" : "test_4", "test_4" : "test_3", "test_5" : "test_5"}
scenarios = []
import glob
gff_files = [f for f in glob.glob(os.path.join(gff_directory, "*.gff")) if os.path.basename(f).replace(".gff", "") in scenario_mapping]
gff_files = list(sorted(gff_files, key=lambda x: scenario_mapping[os.path.basename(x).replace(".gff", "")]))
for gff_file in gff_files:
    file_path = os.path.join(gff_directory, gff_file)
    genes_dict, positions_dict, amr = parse_gff(file_path)
    if os.path.basename(gff_file) != "test_3.gff" and os.path.basename(gff_file) != "test_5.gff":
        genes_dicts.append(genes_dict)
        positions_dicts.append(positions_dict)
    if os.path.basename(gff_file) == "test_3.gff":
        reversed_by_keys = dict(list(reversed(list(genes_dict.items()))))
        genes_dicts.append(reversed_by_keys)
        positions_dicts.append(positions_dict)
        #reversed_by_keys["NC_023289.2:103539-140060_plasmid"] = ["+" + g[1:] if g[0] == "-" else "-" + g[1:] for g in list(reversed(reversed_by_keys["NC_023289.2:103539-140060_plasmid"]))]
    if os.path.basename(gff_file) == "test_5.gff":
        reversed_by_keys = dict(list(reversed(list(genes_dict.items()))))
        genes_dicts.append(reversed_by_keys)
        positions_dicts.append(positions_dict)
    amr_alleles.update(amr)
    scenarios.append(scenario_mapping[os.path.basename(gff_file).replace(".gff", "")])
# Generate unique colors for each gene
unique_genes = set(g[1:] for genes_dict in genes_dicts for genes in genes_dict.values() for g in genes)
gene_colors = {gene: "lightblue" if gene not in amr_alleles else "red" for gene in unique_genes}

out_data = {}
for i in range(len(genes_dicts)):
    scenario = i + 1
    out_data[scenario] = (genes_dicts[i], positions_dicts[i], gene_colors)
with open("aggregated_sim_results/context_plot_inputs.json", "w") as o:
    o.write(json.dumps(out_data))
# Plot the data in batches
plot_reads_gene_positions(genes_dicts, positions_dicts, output_file_base, scenarios, gene_colors)

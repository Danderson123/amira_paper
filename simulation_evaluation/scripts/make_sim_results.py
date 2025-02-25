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

from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from palettable import cartocolors


# Define constants and shared settings
plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 30})

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
        "1": [" "],
        "2": ["chromosome"],
        "3": ["plasmid", "chromosome"],
        "4": ["chromosome"],
        "5": ["plasmid", "chromosome"],
        "6": ["plasmid", "chromosome", "chromosome"]
    }
    if scenario in custom_y_labels:
        ax.set_yticks(range(1, len(custom_y_labels[scenario]) + 1))
        ax.set_yticklabels(custom_y_labels[scenario])
        if scenario != "6":
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

# Main combined plot function
def plot_combined(genes_data, recall_data, output_file):
    depths = [5, 10, 20, 40, 80]
    lengths = [5, 10, 20, 40]
    scenarios = list(genes_data.keys())
    fig = plt.figure(figsize=(30, 45))
    grid = gridspec.GridSpec(5, 2, hspace=0.1, wspace=0.2, width_ratios=[1.5, 1])
    global_max_end = 0
    for s in genes_data:
        _, positions_dict, _ = genes_data[s]
        if len(positions_dict) != 0:
            for contig in positions_dict:
                max_end = max([end for start, end in positions_dict[contig]])
                if max_end > global_max_end:
                    global_max_end = max_end
    for i, scenario in enumerate(scenarios):
        if scenario == "1":
            continue
        genes_dict, positions_dict, gene_colors = genes_data[scenario]
        recalls = recall_data[scenario]
        # Context plot (left column)
        ax1 = fig.add_subplot(grid[i-1, 0])
        plot_gene_positions(ax1, genes_dict, positions_dict, gene_colors, scenario, global_max_end)
        if scenario == "6":
            ax1.set_xlabel("Nucleotide position (bp)")
        # AMR recall plot (right column)
        ax2 = fig.add_subplot(grid[i-1, 1])
        plot_amr_recall(ax2, depths, recalls, lengths, scenario)
        if scenario == "4":
            ax2.set_ylabel("AMR gene recall")
        if scenario == "6":
            ax2.set_xlabel("Read depth (x)")
    # Create a single unified legend for the entire figure
    line_styles = ['-', '--', '-.', ':']
    marker_styles = ['o', 's', '^', 'D']
    legend_colors = generate_colorblind_safe_colors(len(lengths))  # Use the colorblind-safe palette

    legend_elements = [
        plt.Line2D(
            [0], [0],
            linestyle=line_styles[i],
            marker=marker_styles[i],
            color=legend_colors[i],  # Use colors from the palette
            markersize=10,
            label=f"{lengths[i]}kb"
        )
        for i in range(len(lengths))
    ]

    fig.legend(
        handles=legend_elements, 
        loc='center right', 
        bbox_to_anchor=(1.02, 0.5), 
        fontsize=30,
        title="Read Lengths"
    )
    plt.savefig(output_file, dpi=600, bbox_inches='tight')
    plt.savefig(output_file.replace(".png", ".pdf"), bbox_inches='tight')

def adjust_lightness(color, amount=0.5):
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])

def plot_combined_new(genes_data, amira_recall_data, flye_recall_data, output_file):
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
        if scenario == "1":
            continue
        genes_dict, positions_dict, gene_colors = genes_data[scenario]
        recalls_amira = amira_recall_data[scenario]
        recalls_flye = flye_recall_data[scenario]
        # Context plot (left column)
        ax1 = fig.add_subplot(grid[i-1, 0])
        plot_gene_positions(ax1, genes_dict, positions_dict, gene_colors, scenario, global_max_end)
        if scenario == "6":
            ax1.set_xlabel("Nucleotide position (bp)")

        # Lollipop plot for recall differences (right column)
        ax2 = fig.add_subplot(grid[i-1, 1])
        plot_dict = {"Lengths": [], "Depths": [], "Amira": [], "Flye": [], "x": [], "color": []}
        for depth in recalls_amira:
            print(scenario, depth, recalls_amira[depth])
            for i in range(len(recalls_amira[depth])):
                plot_dict["Lengths"].append(int(lengths[i]))
                plot_dict["Depths"].append(int(depth))
                plot_dict["Amira"].append(recalls_amira[depth][i])
                try:
                    plot_dict["Flye"].append(recalls_flye[depth][i])
                except:
                    print(lengths[i], depth)
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
        if scenario == "6":
            ax2.set_xlabel("Mean read length (bp)")
            ax2.set_xticks(xticks_positions)
            ax2.set_xticklabels(custom_labels)  # Custom labels when scenario == "6"
        else:
            ax2.set_xticks(xticks_positions)
            ax2.set_xticklabels([])
        if scenario == "4":
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

with open("aggregated_sim_results/context_plot_inputs.json") as i:
    genes_data = json.load(i)
with open("aggregated_sim_results/recall_plot_inputs.json") as i:
    amira_recall_data = json.load(i)
with open("aggregated_sim_results/flye_recall_plot_inputs.json") as i:
    flye_recall_data = json.load(i)
plot_combined_new(genes_data, amira_recall_data, flye_recall_data, 'aggregated_sim_results/final_combined_plot.png')

import pyfastaq
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import statistics
import subprocess
from Bio import AlignIO
from tqdm import tqdm
import numpy as np
import glob
import json

def rc_sequence(sequence):
    replacement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(list(reversed([replacement[b] for b in list(sequence)])))


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

def write_assembler_alleles(assembler):
    AMRFP_outputs = glob.glob(f"AMR_finder_plus_results.{assembler}_assemblies/*/AMR_finder_plus_results.tsv")
    assembly_dir = f"{assembler}_assemblies"
    output_dir = f"{assembler}_alleles"
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    for o in AMRFP_outputs:
        sample = os.path.basename(os.path.dirname(o))
        if "flye" in assembler or "unicycler" in assembler:
            assembly = os.path.join(assembly_dir, sample, "assembly.fasta")
        if "raven" in assembler:
            assembly = os.path.join(assembly_dir, sample + ".fasta")
        # import the assembly
        reader = pyfastaq.sequences.file_reader(assembly)
        seq_dict = {}
        for sequence in reader:
            seq_dict[str(sequence.id).split(" ")[0]] = str(sequence.seq)
        # import the AMRFP output
        result = pd.read_csv(o, sep="\t")
        gene_seqs = {}
        for index, row in result.iterrows():
            if row["Element subtype"] == "POINT":
                continue
            if "PARTIAL" in row["Method"]:
                continue
            if row["Gene symbol"] not in gene_seqs:
                gene_seqs[row["Gene symbol"]] = {}
            if row["Strand"] == "+":
                gene_seqs[row["Gene symbol"]][f"{row['Gene symbol']}_{len(gene_seqs[row['Gene symbol']])}"] = seq_dict[str(row["Contig id"])][int(row["Start"]) - 1: int(row["Stop"])]
            elif row["Strand"] == "-":
                gene_seqs[row["Gene symbol"]][f"{row['Gene symbol']}_{len(gene_seqs[row['Gene symbol']])}"] = rc_sequence(seq_dict[str(row["Contig id"])][int(row["Start"]) - 1: int(row["Stop"])])
        gene_fastas = []
        for g in gene_seqs:
            for a in gene_seqs[g]:
                gene_fastas.append(f">{a}\n{gene_seqs[g][a]}")
        with open(os.path.join(output_dir, sample + ".fasta"), "w") as outfile:
            outfile.write("\n".join(gene_fastas))

def calculate_allele_accuracy_with_mafft(all_seqs, output_dir):
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
    amira_seqs = [aligned for header, aligned in seqs if "_test" in header]
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
                    max_similarity = similarity_matrix[i, j]
                    best_truth_idx = i
                    best_amira_idx = j
        # If a valid pair was found, mark the truth and amira alleles as paired
        if best_truth_idx != -1 and best_amira_idx != -1:
            paired_similarities.append(max_similarity)
            paired_truths.add(best_truth_idx)
            paired_amiras.add(best_amira_idx)
    return paired_similarities

def plot_nucleotide_results_violin(similarities, output_file):
    # Extract data
    r9_4_1_data = similarities["9.4.1"]
    r10_4_1_data = similarities["10.4.1"]
    # Prepare data for a single plot
    combined_data = []
    all_accuracies = []
    for key, data in zip(["9.4.1", "10.4.1"], [r9_4_1_data, r10_4_1_data]):
        accuracies = []
        for value in data:
            combined_data.append({"Technology": key, "Allele Accuracy": value})
            accuracies.append(value)
            all_accuracies.append(value)
        print(f"Flye {key} accuracy = {statistics.mean(accuracies)}")
    print(f"Flye all allele accuracies = {statistics.mean(all_accuracies)}")
    # Convert to DataFrame for seaborn compatibility
    df = pd.DataFrame(combined_data)
    # Create the plot
    fig = plt.figure(figsize=(12, 12), dpi=600)  # Explicitly set DPI here
    # Violin plot
    sns.violinplot(
        x="Technology",
        y="Allele Accuracy",
        data=df,
        palette="pastel",
        linewidth=1.5,
        saturation=1,
        cut=0
    )
    # Strip plot
    sns.stripplot(
        x="Technology",
        y="Allele Accuracy",
        data=df,
        palette="pastel",
        jitter=True,
        marker='o',
        edgecolor='black',
        linewidth=1,
        alpha=1,
        size=6
    )
    # Customize the plot
    plt.ylabel("Allele accuracy", fontsize=16, labelpad=20)
    plt.ylim([0.6, 1.005])
    plt.xlabel("", fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.grid(axis="y", linestyle="--", alpha=0.7, zorder=1)
    plt.grid(axis="x", visible=False)
    plt.gca().spines['left'].set_visible(True)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['bottom'].set_linewidth(0.5)
    plt.tight_layout()
    # Save the plot
    plt.savefig(output_file, dpi=900)
    plt.savefig(output_file.replace(".png", ".pdf"))
    plt.close()

def calculate_assembler_accuracies(assembler):
    write_assembler_alleles(assembler)
    output_dir = "truth_results"
    # make the output dir
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    # get the sample names
    samples = [os.path.basename(f).replace(".json", "") for f in glob.glob("truth_jsons/*.json")]
    # Initialize dictionary to store data per scenario, depth, and length
    data_list = {s: {} for s in samples}
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
        # import the truth fasta
        with open(os.path.join("truth_allele_sequences", f"{s}.fasta")) as i:
            allele_seqs = i.read().split(">")[1:]
        true_nucleotide_sequences = {}
        for allele in allele_seqs:
            allele_name = allele.split("\n")[0]
            truth_sequence = "".join(allele.split("\n")[1:])
            amira_allele, reference_allele, cellular_copy_number = allele_name.split(";")
            gene_name = apply_rules(reference_allele.split(".")[0])
            if gene_name not in true_nucleotide_sequences:
                true_nucleotide_sequences[gene_name] = []
            true_nucleotide_sequences[gene_name].append(f">{gene_name}_truth\n{truth_sequence}")
        # load the samplea allele fasta
        amrfp_alleles = os.path.join(f"{assembler}_alleles", s + ".fasta")
        # get the amira allele sequences
        test_nucleotide_sequences = {}
        with open(amrfp_alleles) as i:
            test_alleles = i.read().split(">")[1:]
        for row in test_alleles:
            allele_name = row.split("\n")[0]
            gene_name = apply_rules("_".join(allele_name.split("_")[:-1]))
            allele_sequence = "".join(row.split("\n")[1:])
            if gene_name not in test_nucleotide_sequences:
                test_nucleotide_sequences[gene_name] = []
            test_nucleotide_sequences[gene_name].append(f">{gene_name}_test\n{allele_sequence}")
        # get the nucleotide accuracy of the amira alleles
        for gene_name in true_nucleotide_sequences:
            if gene_name in test_nucleotide_sequences:
                all_sequences = "\n".join(true_nucleotide_sequences[gene_name] + test_nucleotide_sequences[gene_name])
                test_similarities = calculate_allele_accuracy_with_mafft(all_sequences, output_dir)
                all_similarities[method] += test_similarities

    return all_similarities
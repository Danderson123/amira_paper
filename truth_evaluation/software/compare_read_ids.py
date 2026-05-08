import gzip
import glob
import os
from tqdm import tqdm

def parse_fastq_lines(fh):
    # Initialize a counter to keep track of the current line number
    line_number = 0
    # Iterate over the lines in the file
    for line in fh:
        # Increment the line number
        line_number += 1
        # If the line number is divisible by 4, it's a sequence identifier line
        if line_number % 4 == 1:
            # Extract the identifier from the line
            identifier = line[1:]
        # If the line number is divisible by 4, it's a sequence line
        elif line_number % 4 == 2:
            sequence = line.strip()
        elif line_number % 4 == 0:
            # Yield the identifier, sequence and quality
            yield identifier, sequence, line.strip()


def parse_fastq(fastq_file):
    # Initialize an empty dictionary to store the results
    results = {}
    # Open the fastq file
    if ".gz" in fastq_file:
        try:
            with gzip.open(fastq_file, "rt") as fh:
                # Iterate over the lines in the file
                for identifier, sequence, quality in parse_fastq_lines(fh):
                    # Add the identifier and sequence to the results dictionary
                    results[identifier.replace("\n", "")] = {
                        "sequence": sequence,
                        "quality": quality,
                    }
            return results
        except OSError:
            pass
    with open(fastq_file, "r") as fh:
        # Iterate over the lines in the file
        for identifier, sequence, quality in parse_fastq_lines(fh):
            # Add the identifier and sequence to the results dictionary
            results[identifier.replace("\n", "")] = {"sequence": sequence, "quality": quality}
    # Return the dictionary of results
    return results

mapping = {
    "AUSMDU00010405": "SRR32405434_1",
    "AUSMDU00015264": "SRR32405436_1",
    "AUSMDU00021208": "SRR32405439_1",
    "AUSMDU00031899": "SRR32405433_1",
    "AUSMDU00031978": "SRR32405440_1",
    "AUSMDU00032793": "SRR32405441_1",
    "AUSMDU00036400": "SRR32405435_1",
    "AUSMDU00040126": "SRR32405437_1",
    "AUSMDU00055259": "SRR32405442_1",
    "AUSMDU00062512": "SRR32405438_1"
}

old_reads = glob.glob("/hps/nobackup/iqbal/dander/Ryan_Wick_E_coli/old_ont_reads/*")
new_reads = glob.glob("/hps/nobackup/iqbal/dander/amira_paper_latest/truth_evaluation/amira_paper_data/truth_evaluation/nanopore_reads/*")

for r1 in tqdm(old_reads):
    old = set([k.split(" ")[0] for k in parse_fastq(r1).keys()])
    r2 = [n for n in new_reads if mapping[os.path.basename(r1).split(".fastq")[0]] in n][0]
    new = set([k.split(" ")[1].split("/")[0] for k in parse_fastq(r2).keys()])
    print(r1, r2, len(new.intersection(old)), len(new), len(old))

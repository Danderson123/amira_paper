from tqdm import tqdm
import subprocess
import os
import requests

# Directories and tracking variables
failed = []
success = []

if not os.path.exists("nanopore_reads"):
    os.mkdir("nanopore_reads")

with open("E_coli_nanopore_accessions.txt") as i:
    nanopore_accessions = i.read().split("\n")

# Filter for nanopore-related sequencing runs
#for index, row in data.iterrows():
#    description = row["description"].lower()
#    if any(s in description for s in {"promethion", "minion", "nanopore", "ont", "pacbio"}):
#        nanopore_accessions.append(row["run_accession"])

# Function to search ENA and retrieve fastq file link
def get_fastq_links(accession):
    url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={accession}&result=read_run&fields=fastq_ftp"
    response = requests.get(url)
    if response.status_code == 200:
        lines = response.text.strip().split("\n")
        if len(lines) > 1:  # Skip header
            fastq_links = lines[1].split("\t")[-1]
            return fastq_links.split(";")  # Return list of fastq URLs
    return None

# Process each accession and download the corresponding fastq files
for accession in tqdm(nanopore_accessions):
    target_file = os.path.join("nanopore_reads", f"{accession}_1.fastq.gz")

    if os.path.exists(target_file):
        print(f"Skipping {accession}, file already exists.")
        success.append(accession)
        continue

    # Search for the fastq links using ENA's API
    fastq_links = get_fastq_links(accession)

    if fastq_links:
        try:
            # Download the first fastq file (assuming paired-end, this can be adjusted)
            fastq_url = f"ftp://{fastq_links[0]}"
            command = f"wget {fastq_url} -O {target_file}"
            subprocess.run(command, shell=True, check=True)
            success.append(accession)
        except Exception as e:
            print(f"Failed to download {accession}: {e}")
            failed.append(accession)
    else:
        print(f"No fastq link found for {accession}.")
        failed.append(accession)

# Write the lists of successful and failed accessions
with open("failed_accessions.txt", "w") as o:
    o.write("\n".join(failed))

with open("successful_accessions.txt", "w") as o:
    o.write("\n".join(success))

import os
import subprocess

def get_mean_read_depth_per_contig(bam_file, samtools_path, core_genes=None):
    # Run samtools depth command and capture output
    command = [samtools_path, "mpileup", "-aa", "-Q", "0", "--ff", "UNMAP,QCFAIL,DUP,SUPPLEMENTARY", bam_file]
    # Run the command and capture the output
    result = subprocess.run(command, capture_output=True, text=True, check=True)
    # Parse the mpileup output
    contig_depths = {}
    contig_positions = {}
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) > 3:  # Ensure sufficient fields exist
            contig = fields[0]  # Contig name is the 1st column
            if core_genes is not None:
                if contig not in core_genes:
                    continue
            depth = int(fields[3])  # 4th column is the depth
            if contig not in contig_depths:
                contig_depths[contig] = 0
                contig_positions[contig] = 0
            contig_depths[contig] += depth
            contig_positions[contig] += 1
    # Calculate mean depth for each contig
    mean_depths = {
        contig: (contig_depths[contig] / contig_positions[contig]) for contig in contig_depths
    }
    print(mean_depths)
    return mean_depths

get_mean_read_depth_per_contig("/hps/nobackup/iqbal/dander/amira_paper_latest/truth_evaluation/amira_output/AUSMDU00032793/AMR_allele_fastqs/longest_reads.bam", "samtools")

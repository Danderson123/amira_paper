import glob
import json
import os

def apply_rules(gene):
    if "blaCTX-M" in gene:
        gene = "blaCTX-M"
    if "blaNDM" in gene:
        gene = "blaNDM"
    if "blaOXA" in gene:
        gene = "blaOXA"
    if "aac6-Ib" in gene:
        gene = "aac6-Ib"
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
    if "aac3-II" in gene and "aac3-III" not in gene:
        gene = "aac3-II"
    if "aac3-III" in gene:
        gene = "aac3-III"
    if "blaSHV" in gene:
        gene = "blaSHV"
    if "qnrS" in gene:
        gene = "qnrS"
    if "aac3-I" in gene and "aac3-II" not in gene and "aac3-III" not in gene:
        gene = "aac3-I"
    if "blaKPC" in gene:
        gene = "blaKPC"
    if "mcr-5" in gene:
        gene = "mcr-5"
    if "qnrB" in gene:
        gene = "qnrB"
    if "cmlA" in gene:
        gene = "cmlA"
    if "aph3-I" in gene and "aph3-II" not in gene and "aph3-III" not in gene:
        gene = "aph3-I"
    if "aph3-II" in gene and "aph3-III" not in gene:
        gene = "aph3-II"
    if "aph3-III" in gene:
        gene = "aph(3)-III"
    if "aph3-I" in gene and "aph3-II" not in gene and "aph3-III" not in gene:
        gene = "aph3-I"
    if "aph3-II" in gene and "aph3-III" not in gene:
        gene = "aph3-II"
    if "aph3-III" in gene:
        gene = "aph3-III"
    if "aph4-I" in gene and "aph4-II" not in gene and "aph4-III" not in gene:
        gene = "aph4-I"
    if "aph4-II" in gene and "aph4-III" not in gene:
        gene = "aph4-II"
    if "aph4-III" in gene:
        gene = "aph4-III"
    if "aph6-I" in gene and "aph6-II" not in gene and "aph6-III" not in gene:
        gene = "aph6-I"
    if "aph6-II" in gene and "aph6-III" not in gene:
        gene = "aph6-II"
    if "aph6-III" in gene:
        gene = "aph6-III"
    if "qacE" in gene:
        gene = "qacE"
    if "blaLAP" in gene:
        gene = "blaLAP"
    if "aac6-I" in gene and "aac6-II" not in gene and "aac6-III" not in gene:
        gene = "aac6-I"
    if "aac6-II" in gene and "aac6-III" not in gene:
        gene = "aac6-II"
    if "aac6-III" in gene:
        gene = "aac6-III"
    if "blaDHA" in gene:
        gene = "blaDHA"
    if "qepA" in gene:
        gene = "qepA"
    if "blaIMI" in gene:
        gene = "blaIMI"
    if "mcr-1" in gene:
        gene = "mcr-1"
    return gene

for j in glob.glob("truth_jsons/*.json"):
    sample = os.path.basename(j).replace(".json","")
    with open(j) as i:
        a = {k.replace(")", "").replace("(", "").replace("'", "") : v for k,v in json.load(i).items()}
    #try:
    with open(f"true_gene_content/{sample}/present_genes.txt") as i:
        r = set([apply_rules(l.split("\t")[0].split(";")[1].split(".NG")[0]) for l in i.read().split("\n")])
    #except:
    #    continue
    print(sample, r - set(a.keys()), set(a.keys()) - r)

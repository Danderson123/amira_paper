import glob
import os
import pandas as pd
import statistics

all_benchmarks = glob.glob("benchmarks/*")
stats = {}
for f in all_benchmarks:
    res = pd.read_csv(f, sep="\t")
    for index, row in res.iterrows():
        runtime = float(row["s"])
        RAM = float(row["max_rss"])
        rule = os.path.basename(f).split(".")[0]
        sample = os.path.basename(f).replace(f"{rule}.", "").replace(".txt", "")
        if rule not in stats:
            stats[rule] = []
        stats[rule].append((sample, runtime, RAM))
for r in stats:
    print(r)
    print(f"Mean runtime: {statistics.mean([t[1] for t in stats[r]])} seconds, minimum: {min([t[1] for t in stats[r]])}, maximum: {max([t[1] for t in stats[r]])}")
    print(f"Mean RAM: {statistics.mean([t[2] for t in stats[r]])} MB, minimum: {min([t[2] for t in stats[r]])}, maximum: {max([t[2] for t in stats[r]])}\n")

# print flye aggregated
flye_stats = stats["flye_assemble"]
amrfp_flye_stats = stats["run_amrfp_on_flye"]
per_sample_runtime = {}
per_sample_RAM = {}
for row in flye_stats:
    per_sample_runtime[row[0]] = row[1]
    per_sample_RAM[row[0]] = row[2]
for row in amrfp_flye_stats:
    per_sample_runtime[row[0]] = row[1]
    per_sample_RAM[row[0]] = max(row[2], per_sample_RAM[row[0]])
print(f"Flye + AMRFP mean runtime: {statistics.mean(list(per_sample_runtime.values()))} seconds, minimum: {min(list(per_sample_runtime.values()))}, maximum: {max(list(per_sample_runtime.values()))}")
print(f"Flye + AMRFP mean RAM: {statistics.mean(list(per_sample_RAM.values()))} MB, minimum: {min(list(per_sample_RAM.values()))}, maximum: {max(list(per_sample_RAM.values()))}\n")

# print shovill aggregated
shovill_stats = stats["shovill_assemble"]
amrfp_shovill_stats = stats["run_amrfp_on_shovill"]
per_sample_runtime = {}
per_sample_RAM = {}
for row in shovill_stats:
    per_sample_runtime[row[0]] = row[1]
    per_sample_RAM[row[0]] = row[2]
for row in amrfp_shovill_stats:
    per_sample_runtime[row[0]] += row[1]
    per_sample_RAM[row[0]] = max(row[2], per_sample_RAM[row[0]])
print(f"shovill + AMRFP mean runtime: {statistics.mean(list(per_sample_runtime.values()))} seconds, minimum: {min(list(per_sample_runtime.values()))}, maximum: {max(list(per_sample_runtime.values()))}")
print(f"shovill + AMRFP mean RAM: {statistics.mean(list(per_sample_RAM.values()))} MB, minimum: {min(list(per_sample_RAM.values()))}, maximum: {max(list(per_sample_RAM.values()))}\n")

# print shovill aggregated
unicycler_stats = stats["unicycler_assemble"]
amrfp_unicycler_stats = stats["run_amrfp_on_unicycler"]
per_sample_runtime = {}
per_sample_RAM = {}
for row in unicycler_stats:
    per_sample_runtime[row[0]] = row[1]
    per_sample_RAM[row[0]] = row[2]
for row in amrfp_unicycler_stats:
    per_sample_runtime[row[0]] += row[1]
    per_sample_RAM[row[0]] = max(row[2], per_sample_RAM[row[0]])
print(f"unicycler + AMRFP mean runtime: {statistics.mean(list(per_sample_runtime.values()))} seconds, minimum: {min(list(per_sample_runtime.values()))}, maximum: {max(list(per_sample_runtime.values()))}")
print(f"unicycler + AMRFP mean RAM: {statistics.mean(list(per_sample_RAM.values()))} MB, minimum: {min(list(per_sample_RAM.values()))}, maximum: {max(list(per_sample_RAM.values()))}\n")
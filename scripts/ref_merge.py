import pandas as pd
import yaml
import os, sys

def info_update(info, new_info):
    if len(info) == 0:
        info.update(new_info)
        return info
    for k in ["library", "max_distance", "min_dist_diff", "barcode_length", "allowed_lowq", "lowq"]:
        assert info[k] == new_info[k], f"Conflict in {k}: {info[k]} vs {new_info[k]}"
    for fate in ["processed", "valid barcode", "invalid barcode (unmatched)", "invalid barcode (ambiguous)", "low quality barcode", "unparseable barcode", "containing N", "recovered"]:
        info["stat"][fate] += new_info["stat"][fate]
    info["stat"]["valid ratio"] = info["stat"]["valid barcode"] / info["stat"]["processed"]

    return info

output_counts = snakemake.output[0]
output_info = snakemake.output[1]

input_counts = []
input_info = []

sys.stdout.write(f"Reading input files: {snakemake.input}\n")

for file in snakemake.input:
    if file.endswith("counts.csv"):
        input_counts.append(file)
    elif file.endswith("barcode_count_info.yaml"):
        input_info.append(file)
    else:
        raise ValueError(f"Unexpected file: {file}")

# input_counts = [os.path.join(snakemake.config["output"], f"barcode_count/{sample}/counts.csv") for sample in samples]
# input_info = [os.path.join(snakemake.config["output"], f"barcode_count/{sample}/barcode_count_info.yaml") for sample in samples]

counts = pd.concat([pd.read_csv(x) for x in input_counts]).groupby("barcode").sum().reset_index()
counts.to_csv(output_counts, index=None)

info = {}
for x in input_info:
    info_update(info, yaml.safe_load(open(x, "r")))

info["stat"]["detected barcodes"] = len(counts)
info["sample"] = snakemake.wildcards["target_ref"]
info["fastq"] = None

yaml.dump(info, open(output_info, "w"))

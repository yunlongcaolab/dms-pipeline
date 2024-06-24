from pathlib import Path
import pandas as pd
import yaml

# .yaml info and .csv scores for all runs in a batch
df_merge = []
info_df = []
sample_files = {}

input_dir = Path(snakemake.config['output']) / "escape_calc" / snakemake.wildcards.batch

for _ in snakemake.input:
    file = Path(_)
    sample = file.parent.name
    if sample not in sample_files:
        sample_files[sample] = {}
    if file.suffix == '.csv':
        sample_files[sample]['scores'] = file
    elif file.suffix == '.yaml':
        sample_files[sample]['info'] = file
    else:
        raise ValueError(f"Unexpected file type: {file}")
    
for sample, files in sample_files.items():
    with open(files['info'], 'r') as f:
        info = yaml.safe_load(f)
    df = pd.read_csv(files['scores']).assign(sample=sample)
    df_merge.append(df)

    info_df.append(pd.DataFrame(info, index=[sample]))

pd.concat(info_df).to_csv(snakemake.output[0],index=False)
pd.concat(df_merge).to_csv(snakemake.output[1],index=False)
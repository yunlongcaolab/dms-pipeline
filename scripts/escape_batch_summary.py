from pathlib import Path
import pandas as pd
import yaml

# .yaml stat and .csv scores for all runs in a batch
df_merge = []
stat_df = []
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
        sample_files[sample]['stat'] = file
    else:
        raise ValueError(f"Unexpected file type: {file}")
    
for sample, files in sample_files.items():
    with open(files['stat'], 'r') as f:
        stat = yaml.safe_load(f)
    df = pd.read_csv(files['scores']).assign(
        sample=sample,
        antibody=stat['antibody'],
    )
    df_merge.append(df)

    stat_df.append(pd.DataFrame(stat, index=[sample]))

pd.concat(stat_df).to_csv(snakemake.output[0],index=False)
pd.concat(df_merge).to_csv(snakemake.output[1],index=False)
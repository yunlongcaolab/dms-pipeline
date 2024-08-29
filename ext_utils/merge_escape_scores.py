import sys, os, yaml
from pathlib import Path

import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input_dir", "-i", help="input directory")
parser.add_argument("--output", "-o", help="output directory")
parser.add_argument("--variant", action='store_true', help="use variant scores instead of single scores [default: False]")
parser.add_argument("--include_failed", action='store_true', help="include failed samples [default: False]")


def main(args):
    # output args
    sys.stderr.write("Arguments:\n")
    for arg in vars(args):
        sys.stderr.write(f"{arg}: {getattr(args, arg)}\n")
    sys.stderr.write("\n")

    os.makedirs(args.output, exist_ok=True)

    if args.variant:
        files = Path(args.input_dir).rglob(f"variant_escape_scores.csv")
        use_cols = ['aa_substitutions', 'escape_score']
    else:
        files = Path(args.input_dir).rglob(f"single_mut_escape_scores.csv")
        use_cols = ['wildtype', 'site', 'mutation', 'mut_escape']

    df = []
    stat_df = []
    files = list(files)

    sys.stderr.write(f"Processing {len(files)} files\n")

    for file in files:
        stat = yaml.safe_load(open(file.parent / 'calc_escape_stat.yaml'))
        info = yaml.safe_load(open(file.parent / 'calc_escape_info.yaml'))

        stat.update(info)
        batch_sample = file.parent.name
        
        if not stat['pass_QC'] and not args.include_failed:
            continue

        df.append(
            pd.read_csv(file)[use_cols].assign(
                sample=batch_sample,
                batch=file.parent.parent.name,
                antibody=stat['antibody'],
                library=stat['library'],
                antigen=stat['antigen'],
            )
        )

        stat_df.append(stat)

    if len(df) > 0:
        pd.concat(df).to_csv(Path(args.output) / 'merge_scores.csv.gz', index=None)
        pd.DataFrame(stat_df).to_csv(Path(args.output) / 'escape_stat.csv.gz', index=None)
    else:
        sys.stderr.write("No data found\n")

if __name__ == "__main__":
    main(parser.parse_args())

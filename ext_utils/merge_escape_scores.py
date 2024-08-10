import sys, os, yaml
from pathlib import Path

import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input_dir", "-i", help="input directory")
parser.add_argument("--output", "-o", default=None, help="output file [default: stdout]")
parser.add_argument("--variant", action='store_true', help="use variant scores instead of single scores [default: False]")
parser.add_argument("--include_failed", action='store_true', help="include failed samples [default: False]")


def main(args):
    if args.output is None:
        out = sys.stdout
    else:
        if not Path(args.output).parent.exists():
            os.makedirs(Path(args.output).parent)
        out = open(args.output, 'w')

    if args.variant:
        files = Path(sys.argv[1]).glob(f"escape_calc/*/*/variant_escape_scores.csv")
        use_cols = ['aa_substitutions', 'escape_score']
    else:
        files = Path(sys.argv[1]).glob(f"escape_calc/*/*/single_mut_escape_scores.csv")
        use_cols = ['wildtype', 'site', 'mutation', 'mut_escape', 'single_mut_escape']

    df = []
    for file in files:
        stat = yaml.safe_load(open(file.parent / 'calc_escape_stat.yaml'))
        batch_sample = file.parent.name
        
        if not stat['pass_QC'] and not args.include_failed:
            continue

        df.append(
            pd.read_csv(file)[use_cols].assign(
                sample=batch_sample,
                library=stat['library'],
                antibody=stat['antibody'],
                batch=file.parent.parent.name,
                pass_QC=stat['pass_QC'],
            )
        )

    if len(df) > 0:
        pd.concat(df).to_csv(out, index=None)
    else:
        sys.stderr.write("No data found\n")

if __name__ == "__main__":
    main(parser.parse_args())
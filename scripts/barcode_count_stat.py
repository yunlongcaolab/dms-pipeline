import yaml
import pandas as pd

df = {
    'sample': [],
    'batch': [],
    'library': [],
    'valid barcode': [],
    'invalid barcode (unmatched)': [],
    'invalid barcode (ambiguous)': [],
    'low quality barcode': [],
    'unparseable barcode': [],
    'containing N': [],
    'recovered': [],
    'valid ratio': [],

    'unique valid barcodes': [],
    'unique invalid barcodes': [],
    'unique low quality barcodes': [],
    'detected mutations (single)': [],
    'detected mutations (any)': []
}

for x in snakemake.input:
    with open(x, 'r') as f:
        y = yaml.safe_load(f)
        for k in ['sample', 'batch', 'library']:
            df[k].append(y[k])
        for k in ['valid barcode', 'invalid barcode (unmatched)', 'invalid barcode (ambiguous)', 'low quality barcode', 'unparseable barcode', 'containing N', 'recovered', 'unique valid barcodes', 'unique invalid barcodes', 'unique low quality barcodes', 'detected mutations (single)', 'detected mutations (any)']:
            df[k].append(y['stat'][k])
        df['valid ratio'].append('%.3f'%(y['stat']['valid ratio']))
df = pd.DataFrame(df)
df.to_csv(snakemake.output[0], index=None)
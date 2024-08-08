import pandas as pd
from pathlib import Path
import yaml

target = snakemake.wildcards.target
merge_name = snakemake.wildcards.library
output_stat = Path(snakemake.output[0]).parent / f'output_stat_info_{merge_name}.yaml'

output_stat_info = {}

df = pd.concat([
    pd.read_csv(f) for f in snakemake.input.tables
]).reset_index(drop=True)

wt_seq = open(snakemake.input.ref).read().strip().split('\n')
assert wt_seq[0].startswith('>'), f'{snakemake.input.ref} is not a fasta file.'
geneseq = ''.join(wt_seq[1:])
assert len(geneseq) % 3 == 0, f'Length of reference sequence in {snakemake.input.ref} is not a multiple of 3.'

if df['target'].nunique() > 1:
    raise ValueError('Multiple targets detected in input files')
if df['target'].to_list()[0] != target:
    raise ValueError(f'Target mismatch in input files. Expected {target}, got {df["target"].to_list()[0]}')


# ambiguous_barcodes = df.groupby('barcode').filter(lambda x: len(x) > 1)['barcode'].unique().tolist()
ambiguous_barcodes = set()
df_out = df.sort_values('barcode').reset_index(drop=True)
lst = 0
for i in range(1, len(df_out)):
    if df_out['barcode'][i] == df_out['barcode'][i-1] and str(df_out['codon_substitutions'][i]) != str(df_out['codon_substitutions'][i-1]):
        ambiguous_barcodes.add(df_out['barcode'][i])

output_stat_info['n_ambiguous_barcodes'] = len(ambiguous_barcodes)
print("Removed ambiguous barcodes:", len(ambiguous_barcodes))

table = df.assign(library=merge_name).query('barcode not in @ambiguous_barcodes').groupby('barcode').agg({
    'target': 'first',
    'library': 'first',
    'variant_call_support': 'sum', 
    'codon_substitutions': 'first',
    'aa_substitutions': 'first',
    'n_codon_substitutions': 'max',
    'n_aa_substitutions': 'max'
}).reset_index()

print("Resulting Complexity:", len(df_out))
output_stat_info['valid_barcodes_after_merging'] = len(df_out)

table.to_csv(snakemake.output[0], index=False)

# count detected single mutations
detected_single_muts = set()
wildtype_aa = {}
for muts in table['aa_substitutions']:
    if pd.notnull(muts) and muts != '':
        for mut in muts.split(' '):
            _site = int(mut[1:-1])
            _wt = mut[0]
            detected_single_muts.add(mut)
            wildtype_aa[_site] = _wt

# Missed sites
missed_sites = []
for site in range(1, len(geneseq) // 3 + 1):
    if site not in wildtype_aa:
        missed_sites.append(site)

output_stat_info['n_detected_sites'] = len(wildtype_aa)
output_stat_info['n_missed_sites'] = len(missed_sites)
output_stat_info['missed_sites'] = list(missed_sites)

print(f"Detected {len(table)} valid barcodes")
print(f"Detected {len(detected_single_muts)} single mutations at {len(wildtype_aa)} sites")

missing_muts = []
for site, ref_aa in wildtype_aa.items():
    for mut_aa in 'ACDEFGHIKLMNPQRSTVWY':
        if mut_aa == ref_aa:
            continue
        mut = f'{ref_aa}{site}{mut_aa}'
        if mut not in detected_single_muts:
            missing_muts.append(mut)

print("Missing single mutations:", missing_muts)

output_stat_info['n_missed_single_mutations'] = len(missing_muts)
output_stat_info['missed_single_mutations'] = missing_muts

# ratio of single-mut variants
single_mut_ratio = len(table.query('n_aa_substitutions == 1')) / len(table)
print(f"Fraction of single-mut variants: {single_mut_ratio:.2f}")
output_stat_info['single_mut_ratio'] = single_mut_ratio

# ratio of WT
wt_ratio = len(table.query('n_aa_substitutions == 0')) / len(table)
print(f"Fraction of WT: {wt_ratio:.2f}")
output_stat_info['WT_ratio'] = wt_ratio

# ratio of multi-mut variants
multi_mut_ratio = len(table.query('n_aa_substitutions > 1')) / len(table)
print(f"Fraction of multi-mut variants: {multi_mut_ratio:.2f}")
output_stat_info['multi_mut_ratio'] = multi_mut_ratio

output_dir = Path(snakemake.output[0]).parent
yaml.dump(output_stat_info, open(output_stat, 'w'))
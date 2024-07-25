import pandas as pd

target = snakemake.wildcards.target
merge_name = snakemake.wildcards.library

df = pd.concat([
    pd.read_csv(f) for f in snakemake.input
])

ambiguous_barcodes = set()
df_out = df.sort_values('barcode').reset_index(drop=True)
lst = 0
for i in range(1, len(df_out)):
    if df_out['barcode'][i] == df_out['barcode'][i-1] and str(df_out['codon_substitutions'][i]) != str(df_out['codon_substitutions'][i-1]):
        ambiguous_barcodes.add(df_out['barcode'][i])
df_out = df_out.assign(library=merge_name).query('barcode not in @ambiguous_barcodes').groupby('barcode').agg({
    'target': 'first',
    'library': 'first',
    'variant_call_support': 'sum', 
    'codon_substitutions': 'first',
    'aa_substitutions': 'first',
    'n_codon_substitutions': 'max',
    'n_aa_substitutions': 'max'
}).reset_index()

print("Removed ambiguous barcodes:", len(ambiguous_barcodes))
print("Resulting Complexity:", len(df_out))
detected_single_muts = set()
wildtype_aa = {}
for muts in df_out['aa_substitutions']:
    if pd.notnull(muts) and muts != '':
        for mut in muts.split(' '):
            _site = int(mut[1:-1])
            _wt = mut[0]
            detected_single_muts.add(mut)
            wildtype_aa[_site] = _wt
min_site = min(wildtype_aa)
max_site = max(wildtype_aa)

print(f"Detected single mutations at {len(detected_single_muts)} sites from {min_site} to {max_site}")

missing_muts = []
for site in range(min_site, max_site + 1):
    if site not in wildtype_aa:
        print('Missing site:', site)
        continue
    ref_aa = wildtype_aa[site]
    for mut_aa in 'ACDEFGHIKLMNPQRSTVWY':
        if mut_aa == ref_aa:
            continue
        mut = f'{ref_aa}{site}{mut_aa}'
        if mut not in detected_single_muts:
            missing_muts.append(mut)

print("Missing single mutations:", missing_muts)

df_out.to_csv(snakemake.output[0], index=None)
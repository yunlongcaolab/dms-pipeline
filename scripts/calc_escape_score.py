import pandas as pd
import numpy as np
import sys, os
from pathlib import Path
import time
import Bio.SeqIO
from Bio.Align import PairwiseAligner
import yaml
import dms_variants.globalepistasis
import dms_variants.codonvarianttable
from binarymap.binarymap import BinaryMap
import collections, itertools
from typing import Dict, Optional, TextIO, Literal

def scaling(x, mn, mx):
    if x < mn:
        return 0
    elif x > mx:
        return 1
    else:
        return (x - mn) / (mx - mn)

def filter_expr_bind(
        df: pd.DataFrame, 
        single_filters: Dict[str, Path], 
        variant_filters: Dict[str, Path], 
        min_single: Dict[str, float], 
        min_variant: Dict[str, float],
        keep_missing: bool, 
        filter_on_variant: bool, 
        filter_on_single: bool,
        log_handle: Optional[TextIO] = sys.stdout) -> pd.DataFrame:

    log_handle.write(f"\nFiltering bad mutants... Start: {time.ctime()}\n")
    escape_scores = df.assign(pass_filters = True)

    # try to apply single mut filters
    if filter_on_single and single_filters is not None:
        for prop, thres in min_single.items():
            if prop not in single_filters:
                log_handle.write(f"Missing filter for {prop}. Skip {prop} single mut filter.\n")
                continue
            if not os.path.exists(single_filters[prop]):
                log_handle.write(f"Missing file {single_filters[prop]}. Skip {prop} single mut filter.\n")
                continue

            single_mut_filt = pd.read_csv(single_filters[prop])

            required_cols = ['mut_start', prop]
            assert set(required_cols).issubset(single_mut_filt.columns), f"Missing columns in {single_filters[prop]}: {required_cols}"

            assert single_mut_filt['mut_start'].nunique() == len(single_mut_filt) # check to make sure each mutation is unique

            if keep_missing: # use negative selection
                log_handle.write("Apply mut_bind_expr filter with negative selection (keep mutants with missing values).\n")
                muts_exclude = set(
                    single_mut_filt.query(
                        f"{prop} < {thres}")['mut_start']
                )
                log_handle.write(str(len(muts_exclude))+" of "+str(len(single_mut_filt))+" mutations have inadequate "+prop+"\n")
                escape_scores[f"muts_pass_{prop}_filter"] = (
                    escape_scores['aa_substitutions'].map(lambda s: set(s.split()).isdisjoint(muts_exclude))
                )
            else: # use positive selection
                log_handle.write("Apply mut_bind_expr filter with positive selection (exclude mutants with missing values).\n")
                muts_adequate = set(
                    single_mut_filt.query(
                        f"{prop} >= {min_single['prop']}")['mut_start']
                )
                log_handle.write(str(len(muts_adequate))+" of "+str(len(single_mut_filt))+" mutations have adequate "+prop+"\n")
                escape_scores[f"muts_pass_{prop}_filter"] = (
                    escape_scores['aa_substitutions'].map(lambda s: set(s.split()).issubset(muts_adequate))
                )
            escape_scores['pass_filters'] = escape_scores['pass_filters'] & escape_scores['muts_pass_'+prop+'_filter']
            log_handle.write(f"{sum(escape_scores['muts_pass_'+prop+'_filter'])} of {len(escape_scores)} mutations pass {prop} filter\n")
    else:
        log_handle.write("Skip single mut filter\n")
    
    # try to apply variant filters
    if filter_on_variant and variant_filters is not None:
        for prop, thres in min_variant.items():
            if prop not in variant_filters:
                log_handle.write(f"Missing filter for {prop}. Skip {prop} variant filter.\n")
                continue

            filt_files = variant_filters[prop]
            var_filt_df_merge = []
            for filt_file in filt_files:
                if not os.path.exists(filt_file):
                    log_handle.write(f"Missing file {filt_file}. Skip {prop} variant filter.\n")
                    continue
                var_filt_df = pd.read_csv(filt_file)
                
                required_cols = ['barcode', prop]
                assert set(required_cols).issubset(var_filt_df.columns), f"Missing columns in {filt_file} for prop {prop}: {required_cols}"
                
                for col in ['library', 'target']:
                    if col in var_filt_df.columns:
                        assert var_filt_df[col].nunique() == 1, f"Multiple {col} in variant filter {prop} {filt_file}"

                var_filt_df_merge.append(var_filt_df)
            
            var_filt_df = pd.concat(var_filt_df_merge, ignore_index=True).groupby('barcode')[prop].mean().reset_index() # merge multiple filters

            log_handle.write(f"Apply variant filter {prop}. Keep missing? {keep_missing}\n")
            variant_pass_df = var_filt_df[['barcode', prop]].dropna().assign(
                **{f'variant_pass_{prop}_filter': lambda x: x[prop] >= thres})
            
            escape_scores = escape_scores.merge(variant_pass_df, how='left', on='barcode')
            escape_scores[f'variant_pass_{prop}_filter'] = escape_scores[f'variant_pass_{prop}_filter'].fillna(
                keep_missing)
            escape_scores['pass_filters'] = escape_scores['pass_filters'] & escape_scores[f'variant_pass_{prop}_filter']
            log_handle.write(f"{sum(escape_scores[f'variant_pass_{prop}_filter'])} of {len(escape_scores)} variants pass {prop} filter\n")
    else:
        log_handle.write("Skip variant filters\n")

    log_handle.write(f"{sum(escape_scores['pass_filters'])} of {len(escape_scores)} mutants pass all filters.\n")
    log_handle.write(f"End filtering: {time.ctime()}\n\n")
    return escape_scores

def calc_epistatsis_model(
        variant_escape_scores: pd.DataFrame, 
        label_sites: Dict[int, str], 
        epistasis_model: str,
        single_mut_scores: Dict[str, float],
        log_handle: Optional[TextIO] = sys.stdout, **kwargs) -> pd.DataFrame:
    log_handle.write(f'\nCalculate for epistasis ... Start: {time.ctime()}\n')
    binary_map = BinaryMap(variant_escape_scores, func_score_col='escape_score')
    if epistasis_model == 'dms_variants':
        # from J. Bloom Lab. https://jbloomlab.github.io/dms_variants/dms_variants.globalepistasis.html
        log_handle.write('Using dms_variants model.\n')
        if 'likelihood' in kwargs:
            if kwargs['likelihood'] == 'Gaussian':
                model = dms_variants.globalepistasis.MonotonicSplineEpistasisGaussianLikelihood(binary_map)
            elif kwargs['likelihood'] == 'Cauchy':
                model = dms_variants.globalepistasis.MonotonicSplineEpistasisCauchyLikelihood(binary_map)
            else:
                raise ValueError(f"Unknown likelihood {kwargs['likelihood']} for dms_variants models.")
        else: # default to Gaussian likelihood
            model = dms_variants.globalepistasis.MonotonicSplineEpistasisGaussianLikelihood(binary_map)
    elif epistasis_model == 'dms_variants_linear':
        log_handle.write('Using dms_variants linear (no epistasis) model.\n')
        model = dms_variants.globalepistasis.NoEpistasisGaussianLikelihood(binary_map)
    else:
        raise ValueError(f"Unknown epistasis model {epistasis_model}.")
    
    model.fit()
    log_handle.write('Epistasis model fitting finished.\n')
    variant_counter = collections.Counter(model.binarymap.substitution_variants)
    muts_counter = collections.Counter(itertools.chain.from_iterable(s.split() for s in model.binarymap.substitution_variants))
    
    effects_df = pd.DataFrame({'mutation': model.binarymap.all_subs}).assign(
        n_single_mut_measurements=lambda x: x['mutation'].map(variant_counter),
        n_any_mut_measurements=lambda x: x['mutation'].map(muts_counter),
    ).pipe(model.add_phenotypes_to_df,  
           substitutions_col='mutation').drop(columns='latent_phenotype').rename(
               columns={'observed_phenotype': 'epistasis_model_score'})
    
    log_handle.write("Removing mutations that do not have EITHER >= 1 single-mutant measurements or >= 3 any-mutant measurements.\n")
    effects_df = effects_df.assign(
        sufficient_measurements=lambda x: (
                    (x['n_single_mut_measurements'] >= 1) |
                    (x['n_any_mut_measurements'] >= 3))
                ).query('sufficient_measurements == True').drop(columns='sufficient_measurements')

    # scaling again
    def score_scaling(x, lower_quantile, upper_quantile):
        mn = max(0,np.nanquantile(x, lower_quantile))
        mx = min(1,np.nanquantile(x, upper_quantile))
        if mx < 10 * mn:
            mx = max(x)
            if mx < 10 * mn:
                mx = 1
        return [scaling(i, mn, mx) for i in x]

    effects_df = effects_df.assign(
        single_mut_escape = lambda x: x['mutation'].map(single_mut_scores),
        site=lambda x: x['mutation'].str[1: -1].astype(int),
        wildtype=lambda x: x['mutation'].str[0]).assign(mutation=lambda x: x['mutation'].str[-1]).assign(
            epistasis_model_score_scaled=lambda x: score_scaling(x['epistasis_model_score'], 0.01, 0.996)).assign(
                # use single if available, otherwise use model score
                mut_escape=lambda x: x['single_mut_escape'].fillna(x['epistasis_model_score_scaled'])
            )

    effects_df['site_number'] = effects_df['site']
    effects_df['site'] = [label_sites[i-1] for i in effects_df['site']]

    site_effects_df = (
        effects_df
        .groupby(['site_number','site'])
        .aggregate(
            site_avg_mut_escape=pd.NamedAgg('mut_escape','mean'),
            site_total_mut_escape=pd.NamedAgg('mut_escape', 'sum'),
            site_avg_single_mut_escape=pd.NamedAgg('single_mut_escape', 'mean'),
            site_total_single_mut_escape=pd.NamedAgg('single_mut_escape', 'sum'),
        )
        .reset_index()
    )

    effects_df = effects_df.merge(site_effects_df, on=['site_number', 'site'], how='left')

    log_handle.write(f'End calculation: {time.ctime()}\n\n')
    return effects_df, site_effects_df

def generate_sequence_numbering(
        ref_seq: str, 
        wt_seq: str, 
        mode: Literal['global', 'local']) -> Dict[int, str]:
    aligner = PairwiseAligner()
    aligner.mode = mode
    aligner.match_score = 3
    aligner.mismatch_score = -1
    aligner.open_gap_score = -4
    aligner.extend_gap_score = -0.5
    aln = aligner.align(ref_seq, wt_seq)[0]
    label_sites = {}
    cur = aln.indices[0, 0]
    offset = 0
    for i in range(aln.shape[1]):
        if aln.indices[0, i] == -1:
            offset += 1
            label_sites[aln.indices[1, i]] = str(cur+1)+'+'+str(offset)
        elif aln.indices[1, i] == -1:
            offset = 0
            cur = aln.indices[0, i]
        else:
            offset = 0
            cur = aln.indices[0, i]
            label_sites[aln.indices[1, i]] = str(cur+1)
    return label_sites

def dms_quality_filter(stat_dict, QC_thres: Dict[str, float], log_handle: Optional[TextIO] = sys.stdout) -> bool:
    log_handle.write(f"\nStart QC ...\n")
    for prop in QC_thres:
        if prop.startswith('min_'):
            prop_name = prop[4:]
            if prop_name in stat_dict:
                if stat_dict[prop_name] < QC_thres[prop]:
                    return False
            else:
                log_handle.write(f"Missing property {prop_name} in stat_dict for QC target {prop}. Skip.\n")

        elif prop.startswith('max_'):
            prop_name = prop[4:]
            if prop_name in stat_dict:
                if stat_dict[prop_name] > QC_thres[prop]:
                    return False
            else:
                log_handle.write(f"Missing property {prop_name} in stat_dict for QC target {prop}. Skip.\n")
        else:
            raise ValueError(f"QC target {prop} must start with 'min_' or 'max_'.")

    return True

# SCRIPT ENTRY
config_output = Path(snakemake.config['output'])

_info = {
    'library': snakemake.params.library,
    'antibody': snakemake.params.antibody,
    'ref': snakemake.params.ref,
    'wt_seq': os.path.abspath(
        snakemake.config['libinfo'][snakemake.params.library]['wt_seq']),
    'ref_numbering_seq': os.path.abspath(
        snakemake.config['libinfo'][snakemake.params.library]['ref_numbering_seq']) if 'ref_numbering_seq' in snakemake.config['libinfo'][snakemake.params.library] else None,
    'table': os.path.abspath(os.path.join(snakemake.config['output'], 'library_tables', 
                          snakemake.config['libinfo'][snakemake.params.library]['target'], 
                          f'{snakemake.params.library}_variant_table.csv')),
    'batch': snakemake.wildcards.batch,
    'sample': snakemake.wildcards.sample,

    'min_ref_count': snakemake.config["calc_escape_scores"]["min_ref_count"],
    'min_variant_support': snakemake.config["calc_escape_scores"]["min_variant_support"],
    'filter_on_single': snakemake.config["calc_escape_scores"]["filter_on_single"],
    'filter_on_variant': snakemake.config["calc_escape_scores"]["filter_on_variant"],
    'include_stop': snakemake.config["calc_escape_scores"]["include_stop"],
    'norm_lower': snakemake.config["calc_escape_scores"]["norm_lower"],
    'norm_upper': snakemake.config["calc_escape_scores"]["norm_upper"],
    # 'normalize_to_WT': snakemake.config["calc_escape_scores"]["normalize_to_WT"],
    'keep_missing': snakemake.config["calc_escape_scores"]["keep_missing"],
    'epistasis': snakemake.config["calc_escape_scores"]["epistasis"],
}

_stat = {}

for info_item in ['library', 'antibody', 'sample']:
    _stat[info_item] = _info[info_item]

for filter_type in ['single', 'variant']:
    if _info[f'filter_on_{filter_type}'] and (snakemake.config["libinfo"][snakemake.params.library][f"{filter_type}_filters"] is not None):
        for prop, files in snakemake.config["libinfo"][snakemake.params.library][f"{filter_type}_filters"].items():
            _info[f'{filter_type}_{prop}_filter'] = [os.path.abspath(x) for x in files] if isinstance(files, list) else os.path.abspath(files)
            try:
                _info[f'{filter_type}_{prop}_min'] = snakemake.config["libinfo"][snakemake.params.library][f"min_{filter_type}"][prop]
            except: # use default value
                _info[f'{filter_type}_{prop}_min'] = snakemake.config["calc_escape_scores"][f"min_{filter_type}"][prop]

table, batch, sample, ref = _info['table'], _info['batch'], _info['sample'], _info['ref']

sample_count = pd.read_csv(config_output / "barcode_count" / sample / "counts.csv")
ref_count = pd.read_csv(config_output / "barcode_count" / ref / "counts.csv")

output_dir = config_output / f"escape_calc/{batch}/{sample}"

# log_handle = open(config_output / snakemake.log.stdout, 'w')
log_handle = sys.stdout

# calculate escape enrichment ratio
ncounts_ref = ref_count['count'].sum().item()
ncounts_escape = sample_count['count'].sum().item()

_stat['ncounts_ref'] = ncounts_ref
_stat['detected_barcodes'] = len(sample_count.query('count > 0'))
_stat['ncounts_escape'] = ncounts_escape

df = ref_count.merge(
    sample_count[['barcode','count']], on = 'barcode', how='left').rename(
        columns = {'count_x':'ref_count','count_y':'sample_count'}
    ).fillna(0).assign(escape_score = lambda x: (x.sample_count/ncounts_escape) / (x.ref_count/ncounts_ref))

min_ref_count = _info['min_ref_count']
min_variant_support = _info['min_variant_support']

df = pd.read_csv(table).drop(columns=['target','library','codon_substitutions','n_codon_substitutions']).merge(df, on = 'barcode', how = 'inner').query(
        'ref_count > @min_ref_count and escape_score >= 0 and variant_call_support >= @min_variant_support'
    )
df['aa_substitutions'] = df['aa_substitutions'].fillna('')
df = df.assign(variant_class = lambda x: x['aa_substitutions'].map(lambda s: 'stop' if '*' in s else 'missense' if len(s.split()) > 0 else 'synonymous'))

_stat['detected_variants'] = len(df.query('sample_count > 0')['aa_substitutions'].unique())
_stat['detected_single_mutants'] = len(df.query('sample_count > 0 and n_aa_substitutions == 1')['aa_substitutions'].unique())

WT_barcodes = df.query('n_aa_substitutions == 0')

WT_enrichment = (WT_barcodes['sample_count'].sum().item()/ncounts_escape) / (WT_barcodes['ref_count'].sum().item()/ncounts_ref)
_stat['WT_enrichment'] = WT_enrichment

df['escape_enrichment_log10_fold_change'] = np.log10(df['escape_score'] / WT_enrichment)
log_handle.write(f"WT enrichment ratio: {WT_enrichment}\n")

# df.to_csv(output_dir / 'variant_effects.csv', index=False)

assert df['escape_score'].max() > 0, "No escape score is greater than 0."

# scaling
mn = np.nanquantile(df['escape_score'], _info["norm_lower"])
mx = np.nanquantile(df['escape_score'], _info["norm_upper"])

assert mn < mx, f"Lower quantile {mn} is not less than upper quantile {mx}"

# adjust error bound if max is 0
_upper_adj = 0
while mx == 0:
    _upper_adj += 1
    mx = np.nanquantile(df['escape_score'], 1.0 - (1.0-_info['norm_upper']) / 2**_upper_adj)

_stat['scaling_lower'] = mn.item()
_stat['scaling_upper'] = mx.item()

log_handle.write(f"\nEnrichment ratio scaling - min: {mn}, max: {mx}, upper adjusted: {1.0 - (1.0-_info['norm_upper']) / 2**_upper_adj}\n")

df = df.assign(raw_escape_score = df['escape_score'])

df['escape_score'] = df['raw_escape_score'].map(lambda x: scaling(x, mn, mx))

df = filter_expr_bind(
    df = df, 
    single_filters = snakemake.config["libinfo"][_info['library']]["single_filters"],
    variant_filters = snakemake.config["libinfo"][_info['library']]["variant_filters"],
    min_single = snakemake.config["calc_escape_scores"]["min_single"],
    min_variant = snakemake.config["calc_escape_scores"]["min_variant"],
    keep_missing = _info["keep_missing"],
    filter_on_variant = _info["filter_on_variant"],
    filter_on_single = _info["filter_on_single"],
    log_handle = log_handle
)

df.to_csv(output_dir / "variant_escape_scores.csv", index=False)

df_filter = df.query('pass_filters == True')

single_mut_scores = df.query('n_aa_substitutions == 1').groupby('aa_substitutions').aggregate(
    total_sample_count = pd.NamedAgg('sample_count', 'sum'),
    total_ref_count = pd.NamedAgg('ref_count', 'sum')
).assign(
    escape_score = lambda x: ((x['total_sample_count']/ncounts_escape) / (x['total_ref_count']/ncounts_ref)).map(lambda y: scaling(y, mn, mx))
)['escape_score'].to_dict()

if not snakemake.config["calc_escape_scores"]["include_stop"]:
    df_filter = df_filter.query('variant_class != "stop"')

if _info['ref_numbering_seq'] is not None:
    ref_numbering_seq = Bio.SeqIO.read(_info['ref_numbering_seq'], 'fasta').seq
    if (len(ref_numbering_seq) % 3 == 0):
        try:
            ref_numbering_seq = ref_numbering_seq.translate()
            log_handle.write("\nReference seq for numbering is translated to a protein sequence.\n")
            
        except:
            log_handle.write("\nReference seq for numbering is considerred as a protein sequence.\n")
    else:
        log_handle.write("\nReference seq for numbering is a protein sequence.\n")
else:
    ref_numbering_seq = None

wt_seq = Bio.SeqIO.read(_info['wt_seq'], 'fasta').seq
if (len(wt_seq) % 3 == 0):
    try:
        wt_seq = wt_seq.translate()
        log_handle.write("\nWT seq is translated to a protein sequence.\n")
    except:
        log_handle.write("\nWT seq is considerred as a protein sequence.\n")
else:
    log_handle.write("\nWT seq is a protein sequence.\n")

if ref_numbering_seq is None:
    if 'seq_offset' in snakemake.config['libinfo'][snakemake.params.library]:
        seq_offset = snakemake.config['libinfo'][snakemake.params.library]['seq_offset']
    else:
        seq_offset = 0
    log_handle.write(f'\nNo alignment for numbering. Only add an offset {seq_offset}\n')
    label_sites = [str(i + seq_offset) for i in range(len(wt_seq))]
else:
    log_handle.write(f'\nAlignment for numbering.\n')
    label_sites = generate_sequence_numbering(ref_numbering_seq, wt_seq, mode='global')


effects_df, site_effects_df = calc_epistatsis_model(
    variant_escape_scores = df_filter,
    label_sites = label_sites, 
    single_mut_scores = single_mut_scores,
    epistasis_model = _info['epistasis'],
    log_handle=log_handle
)

effects_df.to_csv(output_dir / "single_mut_escape_scores.csv", index=False)
site_effects_df.to_csv(output_dir / "site_escape_scores.csv", index=False)

if "QC" in snakemake.config["calc_escape_scores"]:
    _stat["pass_QC"] = dms_quality_filter(_stat, snakemake.config["calc_escape_scores"]["QC"], log_handle=log_handle)
else:
    log_handle.write("\nNo QC threshold is set. Skip QC.\n")
    _stat['pass_QC'] = True

if (_info["filter_on_variant"] or _info["filter_on_single"]) and snakemake.config["calc_escape_scores"]["calc_no_filter"]:
    effects_df, site_effects_df = calc_epistatsis_model(
        variant_escape_scores = df,
        label_sites = label_sites, 
        single_mut_scores = single_mut_scores,
        epistasis_model = _info['epistasis'],
        log_handle=log_handle
    )
    effects_df.to_csv(output_dir / "single_mut_escape_scores_no_filter.csv", index=False)
    site_effects_df.to_csv(output_dir /"site_escape_scores_no_filter.csv", index=False)

yaml.dump(_info, open(output_dir / 'calc_escape_info.yaml', 'w'))
yaml.dump(_stat, open(output_dir / 'calc_escape_stat.yaml', 'w'))

log_handle.close()
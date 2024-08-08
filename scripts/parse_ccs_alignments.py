# Modified from J. Bloom et al.
# https://github.com/jbloomlab/SARS-CoV-2-RBD_DMS/

import alignparse
from alignparse.ccs import Summaries
from alignparse.targets import Targets
from alignparse.constants import CBPALETTE
import alignparse.consensus
from plotnine import theme, element_blank, ggplot, facet_wrap, geom_point, scale_y_log10, scale_x_log10, geom_text, position_dodge, facet_grid, scale_x_continuous
from plotnine import aes, scale_fill_manual, element_text, geom_bar, geom_histogram, geom_vline, ylab, xlab, scale_color_manual,scale_y_continuous
from plotnine.ggplot import save_as_pdf_pages

import pandas as pd
import numpy as np
from pathlib import Path
import yaml

output_dir = Path(snakemake.output[0]).parent
other_cutoff = 0.02  # group as "other" reasons with <= this frac

target = Targets(seqsfile=snakemake.input.ref,
                  feature_parse_specs=snakemake.input.specs)

primary_target = target.target_names[0]
geneseq = target.get_target(primary_target).get_feature('gene').seq

assert len(target.targets) == 1, "Only one target is supported"

pacbio_runs = pd.DataFrame({
    'library': [snakemake.wildcards.library]*len(snakemake.input.reads),
    'name': [Path(x).stem for x in snakemake.input.reads],
    'run': [Path(x).stem for x in snakemake.input.reads],
    'fastq': snakemake.input.reads,
})

ccs_summaries = Summaries(pacbio_runs, report_col = None)
output_stat_info = {}

plots = []
if ccs_summaries.has_zmw_stats():
    p = ccs_summaries.plot_zmw_stats()
    p = p + theme(panel_grid_major_x=element_blank())  # no vertical grid lines
    plots.append(p)
else:
    print('No ZMW stats available.')

for variable in ['length', 'passes', 'accuracy']:
    if ccs_summaries.has_stat(variable):
        p = ccs_summaries.plot_ccs_stats(variable, maxcol=7, bins=25)
        p = p + theme(panel_grid_major_x=element_blank())  # no vertical grid lines
        plots.append(p)
    else:
        print(f"No {variable} statistics available.")

readstats, aligned, filtered = target.parse_alignment(snakemake.input.aln)

processed_ccs = output_dir / 'processed_ccs.csv.gz'

aligned_df = (
    aligned[primary_target].assign(library = snakemake.wildcards.library)
      .drop(columns=['query_clip5', 'query_clip3'])
      .rename(columns={'barcode_sequence': 'barcode'})
)

filtered_df = (
    filtered[primary_target].assign(library = snakemake.wildcards.library)
        .groupby(['library', 'filter_reason'])
        .size()
        .rename('count')
        .reset_index()
        .assign(tot_reason_frac=lambda x: (x.groupby('filter_reason')['count']
                                        .transform('sum')) / x['count'].sum(),
                filter_reason=lambda x: np.where(x['tot_reason_frac'] > other_cutoff,
                                                    x['filter_reason'],
                                                    'other')
                )
    )
readstats = (
    readstats
    .assign(library = snakemake.wildcards.library,
            category_all_targets=lambda x: x['category'].str.split().str[0],
            target=lambda x: x['category'].str.split(None, n=1).str[1],
            valid=lambda x: x['category_all_targets'] == 'aligned')
    )
aligned_df.to_csv(processed_ccs, index=False)
filtered_df.to_csv(output_dir / 'filtered_ccs.csv.gz', index=False)
readstats.to_csv(output_dir / 'readstats.csv', index=False)

output_stat_info['aligned_ccs'] = len(aligned_df)
output_stat_info['filtered_ccs'] = len(filtered_df)

#################
##### plots #####

ncol = 7

# print(readstats)

p = (
    ggplot(readstats
           .groupby(['category_all_targets', 'valid'])
           .aggregate({'count': 'sum'})
           .reset_index(),
           aes('category_all_targets', 'count', fill='valid')) +
    geom_bar(stat='identity') +
    # facet_wrap('~name', ncol=ncol) +
    theme(axis_text_x=element_text(angle=90),
          figure_size=(3,3),
          panel_grid_major_x=element_blank(),
          legend_position='none',
          ) +
    scale_fill_manual(values=CBPALETTE)
    )
plots.append(p)

p = (
    ggplot(readstats
           .groupby(['library', 'category_all_targets', 'valid'])
           .aggregate({'count': 'sum'})
           .reset_index(), 
           aes('category_all_targets', 'count', fill='valid')) +
    geom_bar(stat='identity') +
    facet_wrap('~ library', nrow=1) +
    theme(axis_text_x=element_text(angle=90),
          figure_size=(4, 4),
          panel_grid_major_x=element_blank(),
          legend_position='none',
          ) +
    scale_fill_manual(values=CBPALETTE)
    )
plots.append(p)

p = (
    ggplot(readstats
           .groupby(['target'])
           .aggregate({'count': 'sum'})
           .reset_index(), 
           aes('target', 'count')) +
    geom_point(stat='identity', size=3) +
    theme(axis_text_x=element_text(angle=90),
          figure_size=(0.3 * readstats['target'].nunique(), 2),
          panel_grid_major_x=element_blank(),
          ) +
    scale_y_log10(name='number of reads')
    )
plots.append(p)

p = (
    ggplot(readstats
           .groupby(['target', 'valid'])
           .aggregate({'count': 'sum'})
           .reset_index()
           .assign(total=lambda x: x.groupby('target')['count'].transform('sum'),
                   frac=lambda x: x['count'] / x['total'],
                   ), 
           aes('target', 'frac', fill='valid')) +
    geom_bar(stat='identity') +
    theme(axis_text_x=element_text(angle=90),
          figure_size=(0.5 * readstats['target'].nunique(), 2),
          panel_grid_major_x=element_blank(),
          ) +
    scale_fill_manual(values=CBPALETTE)
    )
plots.append(p)

nreasons = filtered_df['filter_reason'].nunique()

p = (
    ggplot(filtered_df, aes('filter_reason', 'count')) +
    geom_bar(stat='identity') +
    theme(axis_text_x=element_text(angle=90),
          figure_size=(0.3 * nreasons, 2.5),
          panel_grid_major_x=element_blank(),
          )
    )
plots.append(p)
p = (
    ggplot(filtered_df
           .groupby(['library', 'filter_reason'])
           .aggregate({'count': 'sum'})
           .reset_index(),
           aes('filter_reason', 'count')) +
    geom_bar(stat='identity') +
    facet_wrap('~ library', nrow=1) +
    theme(axis_text_x=element_text(angle=90),
          figure_size=(0.3 * nreasons, 2.5),
          panel_grid_major_x=element_blank(),
          )
    )
plots.append(p)

##### plots done #####
save_as_pdf_pages(plots, filename=output_dir / f'{snakemake.wildcards.library}_process_ccs_plots.pdf')
######################

##########################
##### build variants #####
##########################
plots = []

processed_ccs = aligned_df.assign(
    target = primary_target,
    library = snakemake.wildcards.library,
)
nlibs = processed_ccs['library'].nunique()
ntargets = processed_ccs['target'].nunique()

print(f"Read {len(processed_ccs)} CCSs from {nlibs} libraries and {ntargets} targets.")

error_rate_floor = 1e-7
max_error_rate = 1e-4

processed_ccs = (
    processed_ccs
    .assign(barcode_error=lambda x: np.clip(1 - x['barcode_accuracy'],
                                               error_rate_floor, None),
            gene_error=lambda x: np.clip(1 - x['gene_accuracy'],
                                            error_rate_floor, None)
            ).assign(retained=lambda x: ((x['gene_error'] < max_error_rate) &
                                (x['barcode_error'] < max_error_rate)))
    )


#################
##### plots #####
p = (
 ggplot(processed_ccs
        .melt(value_vars=['barcode_error', 'gene_error'],
              var_name='feature_type', value_name='error rate'),
        aes('error rate')) +
 geom_histogram(bins=25) +
 geom_vline(xintercept=max_error_rate,
            linetype='dashed',
            color=CBPALETTE[1]) +
 facet_wrap('~ feature_type') +
 theme(figure_size=(4.5, 2)) +
 ylab('number of CCSs') +
 scale_x_log10()
)
plots.append(p)

p = (
 ggplot(processed_ccs.assign(xlabel=lambda x: x['target'].astype(str) + ', ' + x['library'].astype(str))
                     .groupby(['xlabel', 'retained'])
                     .size()
                     .rename('count')
                     .reset_index(),
        aes('xlabel', 'count', color='retained', label='count')) +
 geom_point(size=3) +
 geom_text(va='bottom', size=7, ha='center',format_string='{:.3g}', nudge_y=0.2) +
 theme(figure_size=(0.5 * nlibs * ntargets, 3),
       panel_grid_major_x=element_blank(),
       axis_text_x=element_text(angle=90),
       ) +
 scale_y_log10(name='number of CCSs') +
 xlab('') +
 scale_color_manual(values=CBPALETTE[1:])
)
plots.append(p)

max_count = 8 # in plot, group all barcodes with >= this many counts
p = (
 ggplot(
    processed_ccs
     .query('retained')
     .groupby(['library', 'barcode'])
     .size()
     .rename('nseqs')
     .reset_index()
     .assign(nseqs=lambda x: np.clip(x['nseqs'], None, max_count)),
    aes('nseqs')) +
 geom_bar() +
 facet_wrap('~ library', nrow=1) +
 theme(figure_size=(1.75 * nlibs, 2),
       panel_grid_major_x=element_blank(),
       ) +
 ylab('number of barcodes') +
 xlab('CCSs for barcode')
)
plots.append(p)

barcode_counts = processed_ccs.assign(barcode_counts=1).groupby(['library', 'barcode'])['barcode_counts'].count().clip(None, max_count).reset_index()

p = (
 ggplot(
    processed_ccs.merge(barcode_counts, on=['library', 'barcode']),
    aes('barcode_counts')) +
 geom_bar() +
 facet_wrap('~ library', nrow=1) +
 theme(figure_size=(1.75 * nlibs, 2),
       panel_grid_major_x=element_blank(),) +
 ylab('number of CCSs') +
 xlab('CCSs per barcode')
)
plots.append(p)
##### plots done #####

processed_ccs = alignparse.consensus.add_mut_info_cols(processed_ccs,
                                                       mutation_col='gene_mutations',
                                                       n_indel_col='n_indels')

processed_ccs = processed_ccs.assign(has_indel=lambda x: x['n_indels'] > 0)
high_acc = max_error_rate / 10
empirical_acc = []

for desc, query_str in [
        ('retained', 'retained'),
        ('retained, no indel', 'retained and not has_indel'),
        ('10X accuracy',
         f"(gene_error < {high_acc}) and (barcode_error < {high_acc})"),
        ('10X accuracy, no indel',
         f"(gene_error < {high_acc}) and (barcode_error < {high_acc}) and not has_indel")
        ]:
    # get just CCSs in that category
    df = processed_ccs.query(query_str)

    if len(df) == 0:
        print(f"No CCSs for {desc}")
        continue

    # compute empirical accuracy
    empirical_acc.append(
        alignparse.consensus.empirical_accuracy(df,
                                                mutation_col='gene_mutations')
        .assign(description=desc)
        .merge(df
               .groupby('library')
               .size()
               .rename('number_CCSs')
               .reset_index()
               )
        )

# make description categorical to preserve order, and annotate as "actual"
# the category ("retained, no indel") that we will use for building variants.
empirical_acc = (
    pd.concat(empirical_acc, ignore_index=True, sort=False)
    .assign(description=lambda x: pd.Categorical(x['description'],
                                                 x['description'].unique(),
                                                 ordered=True),
            actual=lambda x: np.where(x['description'] == 'retained, no indel',
                                         True, False),
            )
    )

consensus, dropped = alignparse.consensus.simple_mutconsensus(
                        processed_ccs.query('retained'),
                        group_cols=('library', 'barcode', 'target'),
                        mutation_col='gene_mutations',
                        )
consensus = alignparse.consensus.add_mut_info_cols(
                    consensus,
                    mutation_col='gene_mutations',
                    sub_str_col='substitutions',
                    n_indel_col='number_of_indels',
                    overwrite_cols=True)


##### plots #####
p = (
 ggplot(processed_ccs,
        aes('retained', fill='has_indel')) +
 geom_bar(position='dodge') +
 geom_text(aes(label='..count..'), stat='count', va='bottom', size=7,
           position=position_dodge(width=0.9), format_string='{:.2g}') +
 theme(figure_size=(2.5 * nlibs, 3),
       panel_grid_major_x=element_blank(),
       ) +
 ylab('number of CCSs') +
 scale_fill_manual(values=CBPALETTE[1:]) +
 facet_wrap('~ library', nrow=1)
 )
plots.append(p)

p = (
    ggplot(empirical_acc,
           aes('description', 'accuracy', color='actual', label='accuracy')
           ) +
    geom_point(size=3) +
    geom_text(va='bottom', size=9, format_string='{:.3g}', nudge_y=0.003) +
    facet_wrap('~ library') +
    theme(figure_size=(1.75 * nlibs, 2.25),
          axis_text_x=element_text(angle=90),
          panel_grid_major_x=element_blank(),
          ) +
    xlab('') +
    scale_y_continuous(name='empirical accuracy', limits=(0.95, 1.005)) +
    scale_color_manual(values=CBPALETTE)
    )
plots.append(p)

max_variant_call_support = 6  # group variants with >= this much support

p = (
 ggplot(consensus
        .assign(variant_call_support=lambda x: np.clip(x['variant_call_support'],
                                                          None,
                                                          max_variant_call_support),
                indel_state=lambda x: np.where(x['number_of_indels'] > 0,
                                                  'has indel', 'no indel')
                ),
        aes('variant_call_support')) +
 geom_bar() +
 ylab('number of variants') +
 facet_grid('indel_state ~ library') +
 theme(figure_size=(1.75 * nlibs, 3.5),
       panel_grid_major_x=element_blank(),
       ) 
 )
plots.append(p)
##### plots done #####

output_stat_info['consensus_barcodes'] = len(consensus)

consensus = consensus.query('number_of_indels == 0')

output_stat_info['consensus_barcodes_remove_indels)'] = len(consensus)

lib_target_counts = (
    consensus
    .groupby(['library', 'target'])
    .size()
    .rename('consensus sequences')
    .reset_index()
    )
p = (ggplot(lib_target_counts.assign(xlabel=lambda x: x['target'] + ', ' + x['library']),
            aes('xlabel', 'consensus sequences')) +
     geom_point(size=3) +
     theme(figure_size=(0.5 * nlibs * ntargets, 2),
           axis_text_x=element_text(angle=90)) +
     xlab('') +
     scale_y_log10()
     )
plots.append(p)

consensus = (
    consensus
    .assign(has_substitutions=lambda x: x['substitutions'].str.len().astype(bool))
    )

has_subs_by_target = (
        consensus
        .groupby(['target', 'library', 'has_substitutions'])
        .aggregate(n_barcodes=pd.NamedAgg('barcode', 'count'))
        .reset_index()
        )

p = (ggplot(has_subs_by_target.assign(xlabel=lambda x: x['target'] + ', ' + x['library']),
            aes('xlabel', 'n_barcodes', color='has_substitutions')) +
     geom_point(size=3, alpha=0.7) +
     theme(figure_size=(0.5 * nlibs * ntargets, 2),
           axis_text_x=element_text(angle=90)) +
     xlab('') +
     scale_y_log10() +
     scale_color_manual(values=CBPALETTE)
     )
plots.append(p)

print(f"Culling the {len(consensus)} barcodes to remove mutated non-primary targets")

consensus = consensus.query('(target == @primary_target) or (has_substitutions == False)')

print(f"Retained {len(consensus)} barcodes after culling")
dup_barcodes = (
    consensus
    .groupby(['library', 'barcode'])
    .size()
    .rename('duplicate_count')
    .reset_index()
    .query('duplicate_count > 1')
    )

print(dup_barcodes.head())
print(f"\nRemoving the {len(dup_barcodes)} duplicated barcodes."
      f"Started with {len(consensus)} barcodes:")

consensus = (
    consensus
    .merge(dup_barcodes, on=['library', 'barcode'], how='outer')
    .query('duplicate_count.isnull()', engine='python')
    )
print(f"After removing duplicates, there are {len(consensus)} barcodes.")

nt_variant_table_file = output_dir / 'consensus.csv.gz'
print(f"Writing nucleotide variants to {nt_variant_table_file}")
consensus.to_csv(nt_variant_table_file, index=False)

max_nseqs = 8  # plot together all barcodes with >= this many sequences

if (len(dropped) > 0):
    p = (
    ggplot(
        dropped.assign(nseqs=lambda x: np.clip(x['nseqs'], None, max_nseqs)),
        aes('nseqs')) + 
    geom_bar() + 
    scale_x_continuous(limits=(1, None)) +
    xlab('number of sequences for barcode') +
    ylab('number of barcodes') +
    facet_grid('library ~ drop_reason') +
    theme(figure_size=(10, 1.5 * nlibs),
        panel_grid_major_x=element_blank(),
        )
    )
    plots.append(p)

print(f"Read gene of {len(geneseq)} nts for {primary_target} from {snakemake.input.ref}")

from dms_variants.codonvarianttable import CodonVariantTable
variants = CodonVariantTable(
                barcode_variant_file=nt_variant_table_file,
                geneseq=geneseq,
                primary_target=primary_target,
                )

# print(variants.n_variants_df(samples=None).pivot_table(index=['target'], columns='library', values='count').head())

max_support = 10  # group variants with >= this much support
p = variants.plotVariantSupportHistogram(max_support=max_support,
                                         widthscale=1.1,
                                         heightscale=0.9)
p = p + theme(panel_grid_major_x=element_blank())  # no vertical grid lines
plots.append(p)
max_muts = 7  # group all variants with >= this many mutations
for mut_type in ['aa', 'codon']:
    p = variants.plotNumMutsHistogram(mut_type, samples=None, max_muts=max_muts,
                                      widthscale=1.1,
                                      heightscale=0.9)
    p = p + theme(panel_grid_major_x=element_blank())  # no vertical grid lines
    plots.append(p)

p = variants.plotNumCodonMutsByType(variant_type='all', samples=None,
                                    ylabel='mutations per variant',
                                    heightscale=0.8)
p = p + theme(panel_grid_major_x=element_blank())  # no vertical grid lines
plots.append(p)
p = variants.plotNumCodonMutsByType(variant_type='all', samples=None,
                                    ylabel='mutations per variant', 
                                    min_support=2, heightscale=0.8)
p = p + theme(panel_grid_major_x=element_blank())  # no vertical grid lines
plots.append(p)

for mut_type in ['aa', 'codon']:
    p = variants.plotMutHeatmap('all', mut_type, samples=None, #libraries='all_only',
                                widthscale=2)
    plots.append(p)


save_as_pdf_pages(plots, filename=output_dir / f'{snakemake.wildcards.library}_build_variants_plots.pdf')
######################

codon_variant_table_file = output_dir / 'variant_table.csv'
print(f"Writing codon-variant table to {codon_variant_table_file}")

table = variants.barcode_variant_df
output_stat_info['valid_barcodes'] = len(table)

table.to_csv(codon_variant_table_file, index=False)

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

output_stat_info['n_detected_single_mutations'] = len(detected_single_muts)

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
with open(output_dir / 'missing_single_mutations.txt', 'w') as f:
    f.write('\n'.join(missing_muts))

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

yaml.dump(output_stat_info, open(output_dir / f'output_stat_info_{snakemake.wildcards.library}.yaml', 'w'))
import altair as alt
import pandas as pd

from pathlib import Path
import sys

def site_plot(data: pd.DataFrame, use_col: str) -> alt.Chart:
    # columns: site_number,site,site_avg_mut_escape,site_total_mut_escape,site_avg_single_mut_escape,site_total_single_mut_escape,sample

    sys.stdout.write(f"Plotting {use_col}...\n")

    base = alt.Chart(data).encode(
            x=alt.X("site_number:N", title="Site"),
            y=alt.X(use_col+":Q", title='site escape'),
            text="site:N",
            tooltip=[
                "site:N", 
                alt.Tooltip(use_col+":Q", format=".3")],
        ).properties(width=800, height=80)

    points = base.mark_circle(color='#66ccff', opacity=0.8)
    text = base.mark_text(align='center')

    return (points + text).facet('sample:N', columns=1).resolve_scale(y='independent').resolve_axis(y='independent')

def stat_plots(stat: pd.DataFrame) -> alt.Chart:
    charts = alt.vconcat()
    # donut chart for column "pass_QC"
    pass_QC = alt.Chart(stat).mark_arc().encode(
        theta='count()',
        color='pass_QC:N'
    ).transform_aggregate(
        count='count()',
        groupby=['pass_QC']
    ).properties(width=200, height=200).facet('library:N', columns=3)
    charts &= pass_QC

    return charts

# SNAKEMAKE_SCRIPT_ENTRY
samples = pd.read_csv(snakemake.input[1])

scores = pd.read_csv(snakemake.input[0]).merge(
    samples[['sample', 'library', 'pass_QC']], on='sample', how='left'
)

pass_df = scores[scores['pass_QC'] == True]
fail_df = scores[scores['pass_QC'] == False]


stat_plots(samples).save(snakemake.output[-1])

for f in snakemake.output[:-1]:
    target = Path(f)
    _status = target.parent.name
    _agg = 'total' if 'total' in target.name else 'avg'
    _model = 'single_mut' if 'single' in target.name else 'mut'
    
    chart = site_plot(
        data=pass_df if _status == 'pass' else fail_df,
        use_col=f"site_{_agg}_{_model}_escape"
    )

    chart.save(f)
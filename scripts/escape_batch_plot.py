import altair as alt
import pandas as pd

def site_plot(data: pd.DataFrame) -> alt.Chart:
    return (
        alt.Chart(data)
        .mark_line(point=True)
        .encode(
            x="site",
            y="score",
            color="sample",
            tooltip=["sample", "site", "score"],
        )
    )

# SNAKEMAKE_SCRIPT_ENTRY
scores = pd.read_csv(snakemake.input[0])
samples = pd.read_csv(snakemake.input[1])


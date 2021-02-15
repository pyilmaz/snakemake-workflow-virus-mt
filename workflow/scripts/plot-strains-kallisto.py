import pandas as pd
import altair as alt

min_fraction = snakemake.params.get("min_fraction", 0.01)

calls = pd.read_csv(snakemake.input[0], sep="\t")

bars = alt.Chart(calls).mark_bar().encode(
    y=alt.Y("target_id:N", title=""),
    x="fraction:Q",
)

text = bars.mark_text(
    align='left',
    baseline='middle',
    dx=3  # Nudges text to right so it doesn't appear on top of the bar
).encode(
    text='est_counts:Q'
)

(bars + text).save(snakemake.output[0])
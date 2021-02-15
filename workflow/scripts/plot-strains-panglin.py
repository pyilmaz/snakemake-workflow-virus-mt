import pandas as pd
import altair as alt

calls = pd.read_csv(snakemake.input[0])
calls = calls.loc[calls["status"] != "fail"]

alt.Chart(calls).mark_bar().encode(
    y=alt.Y("lineage:N", title=""),
    x="probability"
).save(snakemake.output[0])
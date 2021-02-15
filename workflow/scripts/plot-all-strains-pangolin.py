import pandas as pd
import altair as alt

all_calls = []

for f in snakemake.input:
    calls = pd.read_csv(f)
    calls = calls[calls["status"] != "fail"]
    all_calls.append(calls)
    
all_calls = pd.concat(all_calls)
all_calls["count"] = 1

alt.Chart(all_calls).mark_bar().encode(
    x=alt.X("sum(count):Q", title="count"),
    y=alt.Y("lineage", title=""),
    color=alt.Color("probability", scale=alt.Scale(domain=[0.0,1.0]))
).save(snakemake.output[0])
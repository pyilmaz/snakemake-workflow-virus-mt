from collections import defaultdict

import pandas as pd
import altair as alt

mode = snakemake.wildcards.get("mode", "major")

strains = pd.DataFrame(columns=["target_id", "fraction"])

for f in snakemake.input:
    calls = pd.read_csv(f, sep="\t")
    if mode == "major":
        major_strain = calls.iloc[calls["fraction"].idxmax()]
        strains = strains.append(major_strain[["target_id", "fraction"]])
    else:
        strains = strains.append(calls[["target_id", "fraction"]])

strains["count"] = 1
        
alt.Chart(strains).mark_bar().encode(
    x=alt.X("sum(count):Q", axis=alt.Axis(tickMinStep=1)),
    y=alt.Y("target_id:N", title=""),
    color=alt.Color("fraction", bin=True)
).save(snakemake.output[0])
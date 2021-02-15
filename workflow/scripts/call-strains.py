sys.stderr = open(snakemake.log[0], "w")
import pandas as pd
from pathlib import Path

min_fraction = snakemake.params.get("min_fraction", 0.01)

quant = pd.read_csv(Path(snakemake.input[0]) / "abundance.tsv", sep="\t")

# calculate total counts and fraction
total_counts = quant["est_counts"].sum()
quant["fraction"] = quant["est_counts"] / total_counts

# clean up dataframe
quant = quant[["target_id", "fraction", "est_counts"]]

# summarize noise
other = quant.loc[quant["fraction"] < min_fraction, ["fraction", "est_counts"]].sum()
other["target_id"] = "other"
other.name = "other"

# filter dataframe and add noise row
quant = quant.loc[quant["fraction"] >= min_fraction].append(other).set_index("target_id", drop=True)

# store results
quant.to_csv(snakemake.output[0], sep="\t")
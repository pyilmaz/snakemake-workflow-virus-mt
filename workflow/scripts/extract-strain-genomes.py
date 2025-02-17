sys.stderr = open(snakemake.log[0], "w")
genomes_path = snakemake.params.get("save_strains_to", "")

from os import makedirs, path

import pandas as pd


def extract_oldest_sequence_per_strain(
    metadata,
    host_column="Host",
    date_column="Collection date",
    pangolin_lineage_columns="Pango lineage",
    strain_colum="Virus name",
):
    """Get the oldest strain of each pangolin lineage.

    Args:
        metadata (string): Path to the GISAID metadata file

    Returns:
        pd.DataFrame : Dataframe with the oldes strain of each pangolin lineage
    """

    metadata_df = pd.read_csv(metadata, delimiter="\t")

    # select only human hosts
    metadata_df = metadata_df[metadata_df[host_column] == "Human"]

    # select dates with propper format
    metadata_df = metadata_df[metadata_df[date_column].apply(lambda x: len(x)) == 10]

    # selects pangolin lineage that are not None or NaN
    metadata_df = metadata_df[metadata_df[pangolin_lineage_columns] != "None"]
    metadata_df = metadata_df.dropna(subset=[strain_colum, pangolin_lineage_columns])

    # get first instance for occurence of pangolin_lineage
    metadata_df = metadata_df.sort_values(by=[date_column], ascending=True)
    metadata_df.drop_duplicates(
        subset=[pangolin_lineage_columns], keep="first", inplace=True
    )

    return metadata_df[[strain_colum, pangolin_lineage_columns]].sort_values(
        by=[pangolin_lineage_columns]
    )


def get_sequences(strains, sequences, strain_colum="Virus name"):
    """Extract matching sequences from GISAID fasta file

    Args:
        strains (pd.DataFrame) : Dataframe with the oldes strain of each pangolin lineage
        sequences (dict): Path to the GISDAT fasta file

    Returns:
        dict: Dict of matching sequences
    """

    wanted_strains = strains[strain_colum].values
    wanted_sequences = {}
    curr_sequence = None

    with open(sequences, "r") as handle:
        for line in handle.read().splitlines():
            if line.startswith(">"):
                curr_sequence = line
                if curr_sequence.replace(">", "").split("|")[0] in wanted_strains:
                    wanted_sequence = curr_sequence
                    wanted_sequences[wanted_sequence] = ""
                else:
                    wanted_sequence = None
            else:
                if wanted_sequence:
                    wanted_sequences[wanted_sequence] += line

    return wanted_sequences


def write_sequences(
    strain_df,
    wanted_sequences,
    path_to_save_to=genomes_path,
    strain_colum="Virus name",
    pangolin_lineage_columns="Pango lineage",
):
    """Write single genome files and .txt summary for snakemake output

    Args:
        strain_df (pd.DataFrame): Dataframe with the oldes strain of each pangolin lineage
        wanted_sequences (dict): Dict of matching sequences
        path_to_save_to (string, optional): Path where to store the single genomes. Defaults to genomes_path.
    """
    if not path.exists(path_to_save_to):
        makedirs(path_to_save_to)

    look_up_dict = strain_df.set_index(strain_colum).to_dict()[pangolin_lineage_columns]

    strain_genomes = []

    for key in wanted_sequences.keys():
        key, ending = key.replace(">", "").split("|", 1)
        lineage = key.replace(key, look_up_dict[key])
        strain_genomes.append(path_to_save_to + "/" + lineage + ".fasta")
        with open(path_to_save_to + "/" + lineage + ".fasta", "w") as writestream:
            writestream.write(">" + lineage + "\n")
            writestream.write(wanted_sequences[">" + key + "|" + ending] + "\n")

    with open(snakemake.output[0], "w") as snakemake_output:
        for strain in strain_genomes:
            snakemake_output.write("%s\n" % strain)


# get names of strains we want to extract from GISAID fasta file
strain_df = extract_oldest_sequence_per_strain(snakemake.input.metadata)

# extract these strains from the file
wanted_sequences = get_sequences(strain_df, snakemake.input.sequences)

# check if length matches
if len(strain_df) != len(wanted_sequences):
    raise ValueError("Extracted no of sequences do not match.")

# write single strain sequences and txt for snakemake parsing
write_sequences(strain_df, wanted_sequences)

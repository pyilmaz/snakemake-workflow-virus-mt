checkpoint extract_strain_genomes_from_gisaid:
    input:
        metadata=lambda wildcards: config["strain-calling"]["gisaid-metadata"],
        sequences=lambda wildcards: config["strain-calling"]["gisaid-metafasta"],
    output:
        "resources/gisaid/strain-genomes.txt",
    log:
        "logs/extract-strain-genomes.log",
    params:
        save_strains_to=lambda wildcards: config["strain-calling"][
            "extracted-strain-genomes"
        ],
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/extract-strain-genomes.py"


rule cat_genomes:
    input:
        get_strain_genomes,
    output:
        temp("resources/strain-genomes.fasta"),
    log:
        "logs/cat-genomes.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "cat {input} > {output}"


rule kallisto_index:
    input:
        fasta="resources/strain-genomes.fasta",
    output:
        index=temp("resources/strain-genomes.idx"),
    params:
        extra="",
    log:
        "logs/kallisto-index.log",
    threads: 8
    wrapper:
        "0.70.0/bio/kallisto/index"


rule kallisto_quant:
    input:
        fastq=get_reads_after_qc,
        index="resources/strain-genomes.idx",
    output:
        directory("results/{date}/quant/{sample}"),
    params:
        extra="",
    log:
        "logs/{date}/kallisto_quant/{sample}.log",
    threads: 1
    wrapper:
        "0.70.0/bio/kallisto/quant"


rule call_strains_kallisto:
    input:
        quant="results/{date}/quant/{sample}",
        fq1=lambda wildcards: get_reads_after_qc(wildcards, read="1"),
    output:
        "results/{date}/tables/strain-calls/{sample}.strains.kallisto.tsv",
    log:
        "logs/{date}/call-strains/{sample}.log",
    params:
        min_fraction=config["strain-calling"]["min-fraction"],
    conda:
        "../envs/python.yaml"
    notebook:
        "../notebooks/call-strains.py.ipynb"


rule plot_strains_kallisto:
    input:
        "results/{date}/tables/strain-calls/{sample}.strains.kallisto.tsv",
    output:
        report(
            "results/{date}/plots/strain-calls/{sample}.strains.kallisto.svg",
            caption="../report/strain-calls-kallisto.rst",
            category="Kallisto strain calls",
            subcategory="Per sample",
        ),
    log:
        "logs/{date}/plot-strains-kallisto/{sample}.log",
    params:
        min_fraction=config["strain-calling"]["min-fraction"],
    conda:
        "../envs/python.yaml"
    notebook:
        "../notebooks/plot-strains-kallisto.py.ipynb"


rule plot_all_strains_kallisto:
    input:
        lambda wildcards: expand(
            "results/{{date}}/tables/strain-calls/{sample}.strains.kallisto.tsv",
            sample=get_samples_for_date(wildcards.date),
        ),
    output:
        report(
            "results/{date}/plots/all.{mode,(major|any)}-strain.strains.kallisto.svg",
            caption="../report/all-strain-calls-kallisto.rst",
            category="1. Overview",
            subcategory="2. Strain Calls",
        ),
    log:
        "logs/{date}/plot-strains/all.{mode}.log",
    conda:
        "../envs/python.yaml"
    notebook:
        "../notebooks/plot-all-strains-kallisto.py.ipynb"


rule pangolin:
    input:
        contigs="results/{date}/polished-contigs/{sample}.fasta",
        pangoLEARN="results/{date}/pangolin/pangoLEARN",
        lineages="results/{date}/pangolin/lineages",
    output:
        "results/{date}/tables/strain-calls/{sample}.strains.pangolin.csv",
    log:
        "logs/{date}/pangolin/{sample}.log",
    threads: 8
    params:
        pango_data_path=lambda x, input: os.path.dirname(input.pangoLEARN),
    conda:
        "../envs/pangolin.yaml"
    shell:
        "pangolin {input.contigs} --data {params.pango_data_path} --outfile {output} > {log} 2>&1"


rule plot_strains_pangolin:
    input:
        "results/{date}/tables/strain-calls/{sample}.strains.pangolin.csv",
    output:
        report(
            "results/{date}/plots/strain-calls/{sample}.strains.pangolin.svg",
            caption="../report/strain-calls-pangolin.rst",
            category="Pangolin strain calls",
            subcategory="Per sample",
        ),
    log:
        "logs/{date}/plot-strains-pangolin/{sample}.log",
    conda:
        "../envs/python.yaml"
    notebook:
        "../notebooks/plot-strains-pangolin.py.ipynb"


rule plot_all_strains_pangolin:
    input:
        lambda wildcards: expand(
            "results/{{date}}/tables/strain-calls/{sample}.strains.pangolin.csv",
            sample=get_samples_for_date(wildcards.date),
        ),
    output:
        report(
            "results/{date}/plots/all.strains.pangolin.svg",
            caption="../report/all-strain-calls-pangolin.rst",
            category="1. Overview",
            subcategory="2. Strain Calls",
        ),
    log:
        "logs/{date}/plot-strains-pangolin/all.log",
    conda:
        "../envs/python.yaml"
    notebook:
        "../notebooks/plot-all-strains-pangolin.py.ipynb"

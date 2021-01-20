rule freebayes:
    input:
        ref="resources/genome.fasta",
        ref_idx="resources/genome.fasta.fai",
        # you can have a list of samples here
        samples="results/recal/{sample}.bam",
        index="results/recal/{sample}.bam.bai",
    output:
        "results/candidate-calls/{sample}.bcf",
    log:
        "logs/freebayes/{sample}.log",
    params:
        # genotyping is performed by varlociraptor, hence we deactivate it in freebayes by 
        # always setting --pooled-continuous
        extra=(
            "--pooled-continuous --min-alternate-count 1 --min-alternate-fraction 0.01"
        ),
    threads: workflow.cores
    wrapper:
        "0.68.0/bio/freebayes"


rule render_scenario:
    input:
        local("resources/scenario.yaml"),
    output:
        report(
            "results/scenarios/{sample}.yaml",
            caption="../report/scenario.rst",
            category="Variant calling scenarios",
        ),
    log:
        "logs/render-scenario/{sample}.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "sed 's/sample:/{wildcards.sample}:/' {input} > {output}"


rule varlociraptor_preprocess:
    input:
        ref="resources/genome.fasta",
        ref_idx="resources/genome.fasta.fai",
        candidates="results/candidate-calls/{sample}.bcf",
        bam="results/recal/{sample}.bam",
        bai="results/recal/{sample}.bam.bai",
    output:
        "results/observations/{sample}.bcf",
    params:
        depth=config["variant-calling"]["max-read-depth"],
    log:
        "logs/varlociraptor/preprocess/{sample}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor preprocess variants --candidates {input.candidates} "
        "{input.ref} --bam {input.bam} --max-depth {params.depth} --output {output} 2> {log}"


rule varlociraptor_call:
    input:
        obs="results/observations/{sample}.bcf",
        scenario="results/scenarios/{sample}.yaml",
    output:
        "results/calls/{sample}.bcf",
    log:
        "logs/varlociraptor/call/{sample}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor "
        "call variants generic --obs {wildcards.sample}={input.obs} "
        "--scenario {input.scenario} > {output} 2> {log}"

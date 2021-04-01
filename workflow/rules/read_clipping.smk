rule clip_primer:
    input:
        bam=expand(
            "results/{{date}}/mapped/ref~{ref}/{{sample}}.bam",
            ref=config["adapters"]["amplicon-reference"],
        ),
        bed=config["adapters"]["amplicon-primers"],
    output:
        sortbam=temp("results/{date}/clipped-reads/{sample}.bam"),
        sortindex=temp("results/{date}/clipped-reads/{sample}.bam.bai"),
        clippedbam=temp("results/{date}/clipped-reads/{sample}.primerclipped.bam"),
        sortclippedbam=temp(
            "results/{date}/clipped-reads/{sample}.sort.primerclipped.bam"
        ),
        fq1="results/{date}/clipped-reads/{sample}.1.fastq.gz",
        fq2="results/{date}/clipped-reads/{sample}.2.fastq.gz",
    log:
        "logs/{date}/primer-clipping/{sample}.log",
    params:
        dir=lambda w, output: os.path.dirname(output.sortbam),
        bam=lambda w, output: output.sortbam.split("/")[-1],
        dir_depth=lambda w, output: "".join(
            ["../"] * (len(output.sortbam.split("/")) - 1)
        ),
    conda:
        "../envs/bamclipper.yaml"
    threads: 10
    shell:
        """
        samtools sort -@ {threads} -o {output.sortbam} {input.bam} > {log} 2>&1
        samtools index {output.sortbam} >> {log} 2>&1
        cd {params.dir}
        bamclipper.sh -b {params.bam} -p {params.dir_depth}{input.bed} -n {threads} >> {params.dir_depth}{log} 2>&1
        cd {params.dir_depth}
        samtools sort  -@ {threads} -n {output.clippedbam} -o {output.sortclippedbam}  >> {log} 2>&1
        samtools fastq -@ {threads} {output.sortclippedbam} -1 {output.fq1} -2 {output.fq2}  >> {log} 2>&1
        """

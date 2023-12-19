



configfile: "config.yaml"
sample= config["sample"]

rule all:
    input:
        expand("{sample}_bowtie.sorted.bam.bai", sample=config["sample"])




rule bowtie2:
    input:
        R1="{sample}_1.fastq",
        R2="{sample}_2.fastq"
    output:
        "{sample}_bowtie.sam"
    params:
        index=["GRCh38_noalt_as"]
    shell:
        "bowtie2 -x {params.index} -1 {input[0]} -2 {input[1]} -S {output}"



rule samtools_view:
    input:
        "{sample}_bowtie.sam"
    output:
        "{sample}_bowtie.bam"
    shell:
        "samtools view -b -S {input} > {output}"


rule samtools_sort:
    input:
        "{sample}_bowtie.bam"
    output:
        "{sample}_bowtie.sorted.bam"
    shell:
        "samtools sort {input} -o {output}"


rule samtools_index:
    input:
        "{sample}_bowtie.sorted.bam"
    output:
        "{sample}_bowtie.sorted.bam.bai"
    shell:
        "samtools index {input}"

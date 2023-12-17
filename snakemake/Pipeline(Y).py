
    configfile : "/Users/beyzaaygun/Desktop/Bioinformatik/snakemake/yusuf/config.yaml"
rule all:
    input:

rule downloadReference:
    output:
    "Genome_ref/Homo_sapiens.GRCh38.86.gtf"
    "Genome_ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    shell:
        "wget http://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz && \
        wget http://ftp.ensembl.org/pub/release-86/gtf/homo_sapiens/Homo_sapiens.GRCh38.86.gtf.gz && \
        gunzip Homo_sapiens.GRCh38.86.gtf.gz && \
        gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz && \
        mkdir -p Genome_ref && \
        mv Homo_sapiens.GRCh38.86.gtf Genome_ref/ && \
        mv Homo_sapiens.GRCh38.dna.primary_assembly.fa Genome_ref/"



rule get_SRA_by_accession:
"""
Retrieve FASTQ file from SRA (Sequence Read Archive) by run accession number.
"""
    output:
        "{accession}_reads/{accession}_1.fastq",
        "{accession}_reads/{accession}_2.fastq"
    params:
        args = "--split-files --progress --details",
        accession = "{accession}"
    log:
        "{accession}_reads/{accession}.log"
    conda:
        "sra_env.yml"
    shell:
        'mkdir -p {params.accession}_reads && \
        fasterq-dump {params.args} {params.accession} -O {params.accession}_reads'



rule trim:
    input:
        r1="{accession}_reads/{accession}_1.fastq",
        r2="{accession}_reads/{accession}_2.fastq"
    output:
        tr1="RNA-seq/SR1_trimmed.gz",
        tr2="RNA-seq/SR2_trimmed.gz"
    shell:
    "mkdir -p RN-seq/ && \
    fastq-mcf -q 10 {input.r1} {input.r2} -o {output.tr1} -o {output.tr2}"




rule star_index:
    input:
        fasta="Genome_ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        annotation="Genome_ref/Homo_sapiens.GRCh38.86.gtf"
    output:
        "humanINDEX/Genome",
        "humanINDEX/Homo_sapiens.GRCh38.86.gtf",
        "humanINDEX/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        "humanINDEX/SA",
        directory("humanINDEX")
    message:
        "Creating STAR index"
    shell:
        "mkdir -p humanINDEX"
        "STAR --runMode genomeGenerate --genomeDir humanINDEX \
        --genomeFastaFiles {input.fasta} \
        --sjdbGTFfile {input.annotation}"





rule star_align:
    input:
        fq1="{accession}_reads/{accession}_1.fastq",
        fq2="{accession}_reads/{accession}_2.fastq",
        idx="humanINDEX/"
    output:
        '{accession}_STAR/{accession}_Log.final.out'
    params:
        args="--outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 16000000000 --outSAMunmapped Within --twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM"
        outdir = "{accession}_STAR",
        prefix = "{accession}_"
    shell:
        'mkdir -p {params.outdir} && '
        'STAR --genomeDir {input.idx} \
        {params.arg} \
        --readFilesIn {input.fq1},{input.fq2} \
        --outFileNamePrefix {params.outdir}/{params.prefix}'



rule prepare_reference:
    input:
        fasta="Genome_ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        annotation="Genome_ref/Homo_sapiens.GRCh38.86.gtf"
    output:
        seq="rsem_index/reference.seq",
        grp="rsem_index/reference.grp",
        ti="rsem_index/reference.ti"
    params:
        outdir="rsem_index"
    shell:
        "mkdir -p {params.outdir}"
        "rsem-prepare-reference --gtf {input.annotation} {input.fasta} {params.outdir}"
rule ...

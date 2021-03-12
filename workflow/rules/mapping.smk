from multiprocessing import cpu_count

def get_samples_data(wildcards):
    if config["datatype"] == "ILLUMINA":
        paired = (len(config["samples"][wildcards.sample]) == 2)
        assert paired, f"Error in sample {wildcards.sample}, must be paired end."
        return config["samples"][wildcards.sample]
    elif config["datatype"] == "NANOPORE" or config["datatype"] == "PACBIO":
        paired = (len(config["samples"][wildcards.sample]) == 1)
        assert paired, f"Error in sample {wildcards.sample}, must contain a single fasta file."
        return config["samples"][wildcards.sample]


rule minimap2__map_illumina_reads:
    """
    For input preprocessed reads minimap2 finds the most similar genomic region in the provided reference genome.
    Samtools then sort aligned reads according to mapped position on reference genome.
    :input reference: Reference genomic sequences in fasta format
    :input query: Illumina gzipped fastq files with left and right reads e.g. ['path/to/{sample}_R1.fastq.gz','path/to/{sample}_R2.fastq.gz']
    :output bam: Ordered mapped reads according to their location on reference genome
    """
    input:
        reference = config["fasta"],
        query = get_samples_data
    output:
        f"{OUTDIR}/{{sample}}/minimap2/{{sample}}.illumina.sorted.bam"
    log:
        out = f"{OUTDIR}/{{sample}}/logs/minimap2/{{sample}}.illumina.sorted.bam.out",
        err = f"{OUTDIR}/{{sample}}/logs/minimap2/{{sample}}.illumina.sorted.bam.err"    
    message: "Executing {rule}: Aligning {input.query} to {input.reference} using minimap2 and sorting the aligned reads using samtools."
    benchmark:
        f"{OUTDIR}/{{sample}}/logs/minimap2/{{sample}}.illumina.sorted.bam.benchmark"
    threads: 
        lambda cores: cpu_count() - 2
    conda: "../envs/minimap2.yaml"
    shell:
        """
        minimap2 -ax sr -t {threads} \
            -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}" \
            {input.reference} {input.query} | \
        samtools sort -@ {threads} -o {output} - 1> {log.out} 2> {log.err}
        """

rule minimap2__map_nanopore_reads:
    """
    For input preprocessed reads minimap2 finds the most similar genomic region in the provided reference genome.
    Samtools then sort aligned reads according to mapped position on reference genome.
    :input reference: Reference genomic sequences in fasta format
    :input query: Nanopore gzipped fastq file with reads e.g. ['path/to/{sample}.fastq.gz']
    :output bam: Ordered mapped reads according to their location on reference genome
    """
    input:
        reference = config["fasta"],
        query = get_samples_data
    output:
        f"{OUTDIR}/{{sample}}/minimap2/{{sample}}.nanopore.sorted.bam"
    log:
        out = f"{OUTDIR}/{{sample}}/logs/minimap2/{{sample}}.nanopore.sorted.bam.out",
        err = f"{OUTDIR}/{{sample}}/logs/minimap2/{{sample}}.nanopore.sorted.bam.err"
    message: "Executing {rule}: Aligning {input.query} to {input.reference} using minimap2 and sorting the aligned reads using samtools."
    benchmark:
        f"{OUTDIR}/{{sample}}/logs/minimap2/{{sample}}.nanopore.sorted.bam.benchmark"
    threads: 
        lambda cores: cpu_count() - 2
    conda: "../envs/minimap2.yaml"
    shell:
        """
        minimap2 -a -z 600,200 -ax map-ont -t {threads} \
            -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}" \
            {input.reference} {input.query} | \
        samtools sort -@ {threads} -o {output} - 1> {log.out} 2> {log.err}
        """

rule samtools__bam_index:
    """
    Generate .bai index to .bam files to quick recover reads from genomic location of interest.
    :input bam: Mapped reads in bam format
    :output bai: Index of mapped reads to enable fast read retrieval from desired genomic region
    """
    input:
        bam = f"{OUTDIR}/{{sample}}/minimap2/{{sample}}.{{mapper}}.sorted.bam"
    output:
        f"{OUTDIR}/{{sample}}/minimap2/{{sample}}.{{mapper}}.sorted.bam.bai"
    log:
        out = f"{OUTDIR}/{{sample}}/logs/samtools_index/{{sample}}.{{mapper}}.bai.out",
        err = f"{OUTDIR}/{{sample}}/logs/samtools_index/{{sample}}.{{mapper}}.bai.err"
    message: "Executing {rule}: Indexing {input.bam} using samtools."
    benchmark:
        f"{OUTDIR}/{{sample}}/logs/samtools_index/{{sample}}.{{mapper}}.bai.benchmark"
    threads:
        lambda cores: cpu_count() - 2
    conda: "../envs/samtools.yaml"
    shell:
        """
        samtools index -@ {threads} {input.bam} 1> {log.out} 2> {log.err}
        """

rule pbmm2__map_pacbio_reads:
    """
    For input preprocessed reads pbmm2 finds the most similar genomic region in the provided reference genome.
    :input reference: Reference genomic sequences in fasta format
    :input query: PacBio gzipped fastq file with reads e.g. ['path/to/{sample}.fastq.gz']
    :output bam: Ordered mapped reads according to their location on reference genome
    """
    input:
        reference = config["fasta"],
        query     = get_samples_data
    output:
        bam=f"{OUTDIR}/{{sample}}/pbmm2/{{sample}}.pacbio.sorted.bam",
        index=f"{OUTDIR}/{{sample}}/pbmm2/{{sample}}.pacbio.sorted.bai"
    params:
        preset="CCS",
        extra="--sort --unmapped -c 0 -y 70",
        loglevel="INFO"
    log:
        out = f"{OUTDIR}/{{sample}}/logs/pbmm2/{{sample}}.pacbio.sorted.bam.out",
        err = f"{OUTDIR}/{{sample}}/logs/pbmm2/{{sample}}.pacbio.sorted.bam.err"
    message: "Executing {rule}: Aligning {input.query} to {input.reference} using pbmm2 and sorting the aligned reads."
    benchmark:
        f"{OUTDIR}/{{sample}}/logs/pbmm2/{{sample}}.pacbio.sorted.bam.benchmark"
    threads:
        lambda cores: cpu_count() - 2
    conda: "../envs/pbmm2.yaml"
    shell:
        """
        pbmm2 align --num-threads {threads} \
            --preset {params.preset} \
            --rg "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}" \
            --log-level {params.loglevel} \
            {extra} \
            {input.reference} \
            {input.query} \
            {output.bam}) \
        1> {log.out} \
        2> {log.err}
        """

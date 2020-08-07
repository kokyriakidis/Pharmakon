def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = samples.loc[wildcards.sample, ["fq1", "fq2"]].dropna()
    return {"r1": fastqs.fq1, "r2": fastqs.fq2}

samples = pd.read_table(config["samples"]).set_index("sample_name", drop=False)

rule fastqc:
    input:
        unpack(get_fastq)
    output:
        html="qc/fastqc/{sample}/{sample}.html",
        zip="qc/fastqc/{sample}/{sample}.zip"
    wrapper:
        "0.64.0/bio/fastqc"


rule samtools_stats:
    input:
        "dedup/{sample}/{sample}.bam"
    output:
        "qc/samtools-stats/{sample}/{sample}.txt"
    log:
        "logs/samtools-stats/{sample}/{sample}.log"
    wrapper:
        "0.64.0/bio/samtools/stats"


rule multiqc:
    input:
        expand(["qc/samtools-stats/{sample}/{sample}.txt",
                "qc/fastqc/{sample}/{sample}.zip",
                "qc/dedup/{sample}/{sample}.metrics.txt"],
                sample=samples.index)
    output:
        "qc/multiqc/multiqc.html"
    log:
        "logs/multiqc/multiqc.log"
    wrapper:
        "0.64.0/bio/multiqc"

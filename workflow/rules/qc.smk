def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = samples.loc[wildcards.sample, ["fq1", "fq2"]].dropna()
    return {"r1": fastqs.fq1, "r2": fastqs.fq2}

def get_deepvariant_bams(bai=False):
    # case 1: no duplicate removal
    f = "{project_dir}/{sample}/mapped/{sample}.sorted.bam"
    if config["processing"]["remove-duplicates"]:
        # case 2: remove duplicates
        f = "{project_dir}/{sample}/dedup/{sample}.bam"
    return f

def get_gatk_bams(wildcards):
    """Get all aligned reads of given sample."""
    return expand("{project_dir}/{sample}/applybqsr/{sample}.bam",
                  sample=wildcards.sample,
                  project_dir=config["project_dir"])

rule fastqc:
    input:
        unpack(get_fastq)
    output:
        html="{project_dir}/{sample}/qc/fastqc/{sample}.html",
        zip="{project_dir}/{sample}/qc/fastqc/{sample}.zip"
    wrapper:
        "0.64.0/bio/fastqc"

if _platform == "darwin" or config["variant_tool"] == "gatk":

    rule samtools_stats:
        input:
            get_gatk_bams
        output:
            "{project_dir}/{sample}/qc/samtools-stats/{sample}.txt"
        log:
            "{project_dir}/{sample}/logs/samtools-stats/{sample}.log"
        wrapper:
            "0.64.0/bio/samtools/stats"

else:

    rule samtools_stats:
        input:
            get_deepvariant_bams
        output:
            "{project_dir}/{sample}/qc/samtools-stats/{sample}.txt"
        log:
            "{project_dir}/{sample}/logs/samtools-stats/{sample}.log"
        wrapper:
            "0.64.0/bio/samtools/stats"


rule multiqc:
    input:
        expand(["{project_dir}/{sample}/qc/samtools-stats/{sample}.txt",
                "{project_dir}/{sample}/qc/fastqc/{sample}.zip",
                "{project_dir}/{sample}/qc/dedup/{sample}.metrics.txt",
                "{project_dir}/{sample}/recal/{sample}.grp"],
                sample=samples.index,
                project_dir=config["project_dir"])
    output:
        expand("{project_dir}/{sample}/qc/multiqc/multiqc.html", project_dir=config["project_dir"], sample=samples.index)
    log:
        expand("{project_dir}/{sample}/logs/multiqc/multiqc.log", project_dir=config["project_dir"], sample=samples.index)
    wrapper:
        "0.64.0/bio/multiqc"

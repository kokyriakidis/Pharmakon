def get_map_reads_input(wildcards):
    if is_paired_end(wildcards.sample):
        return samples.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()

def is_paired_end(sample):
    sample_unit = samples.loc[sample]
    fq1_null = pd.isnull(sample_unit["fq1"])
    fq2_null = pd.isnull(sample_unit["fq2"])
    paired = not fq1_null and not fq2_null
    assert paired, "Error in sample {}, must be paired end.".format(sample)
    return paired

def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA'".format(sample=wildcards.sample)

if _platform == "darwin":
    rule map_reads:
        input:
            reads=get_map_reads_input,
            idx=rules.get_genome.output
        output:
            "mapped/{sample}/{sample}.sorted.bam"
        log:
            "logs/bwa_mem/{sample}/{sample}.log"
        params:
            index=lambda w, input: input.idx[0],
            extra=get_read_group,
            sort="samtools",
            sort_order="queryname"
        threads: 8
        wrapper:
            "0.64.0/bio/bwa/mem"

elif _platform == "linux" or _platform == "linux2":
    rule map_reads:
        input:
            reads=get_map_reads_input,
            idx=rules.get_genome.output
        output:
            "mapped/{sample}/{sample}.sorted.bam"
        log:
            "logs/bwa_mem2/{sample}/{sample}.log"
        params:
            index=lambda w, input: input.idx[0],
            extra=get_read_group,
            sort="samtools",
            sort_order="queryname"
        threads: 8
        wrapper:
            "0.64.0/bio/bwa-mem2/mem"


rule mark_duplicates:
    input:
        "mapped/{sample}/{sample}.sorted.bam"
    output:
        bam=temp("dedup/{sample}/{sample}.bam"),
        metrics="qc/dedup/{sample}/{sample}.metrics.txt"
    log:
        "logs/picard/dedup/{sample}/{sample}.log"
    params:
        config["params"]["picard"]["MarkDuplicates"]
    wrapper:
        "0.64.0/bio/picard/markduplicates"


rule samtools_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    wrapper:
        "0.64.0/bio/samtools/index"

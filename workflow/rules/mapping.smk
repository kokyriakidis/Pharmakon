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

##### GATK SPECIFIC FUNCTIONS #####

def get_recal_input(bai=False):
    # case 1: no duplicate removal
    f = "{project_dir}/{sample}/mapped/{sample}.sorted.bam"
    if config["processing"]["remove-duplicates"]:
        # case 2: remove duplicates
        f = "{project_dir}/{sample}/dedup/{sample}.bam"
    if bai:
        # case 3: need an index because random access is required
        f += ".bai"
        return f
    else:
        return f

def get_gatk_regions_param(regions=config["processing"].get("restrict-regions"), default=""):
    if regions:
        params = "--intervals '{}' ".format(regions)
        padding = config["processing"].get("region-padding")
        if padding:
            params += "--interval-padding {}".format(padding)
        return params
    return default

def get_recal_params(wildcards, input):
    return (get_gatk_regions_param(regions=input.regions, default="") +
            config["params"]["gatk"]["BaseRecalibrator"])

##### END OF GATK SPECIFIC FUNCTIONS #####


if _platform == "darwin":
    
    rule map_reads:
        input:
            reads=get_map_reads_input
        output:
            "{project_dir}/{sample}/mapped/{sample}.sorted.bam"
        log:
            "{project_dir}/{sample}/logs/bwa_mem/{sample}.log"
        params:
            index="{output_dir}/{genome_build}/bwa/{genome_build}.fa".format(output_dir=config["ref"]["output_dir"], genome_build=config["ref"]["build"]),
            extra=get_read_group,
            sort="samtools",
            sort_order="coordinate"
        threads: 8
        wrapper:
            "0.64.0/bio/bwa/mem"

elif _platform == "linux" or _platform == "linux2":
    
    rule map_reads:
        input:
            reads=get_map_reads_input
        output:
            "{project_dir}/{sample}/mapped/{sample}.sorted.bam"
        log:
            "{project_dir}/{sample}/logs/bwa_mem2/{sample}.log"
        params:
            index="{output_dir}/{genome_build}/bwa/{genome_build}.fa".format(output_dir=config["ref"]["output_dir"], genome_build=config["ref"]["build"]),
            extra=get_read_group,
            sort="samtools",
            sort_order="coordinate"
        threads: 8
        wrapper:
            "0.64.0/bio/bwa-mem2/mem"


rule mark_duplicates:
    input:
        "{project_dir}/{sample}/mapped/{sample}.sorted.bam"
    output:
        bam="{project_dir}/{sample}/dedup/{sample}.bam",
        metrics="{project_dir}/{sample}/qc/dedup/{sample}.metrics.txt"
    log:
        "{project_dir}/{sample}/logs/picard/dedup/{sample}.log"
    params:
        config["params"]["picard"]["MarkDuplicates"]
    wrapper:
        "0.64.0/bio/picard/markduplicates"


rule get_selected_intervals:
    input:
        expand("{project_dir}/{sample}/dedup/{sample}.bam", project_dir=config["project_dir"], sample=samples.index) if config["processing"]["remove-duplicates"] == True else expand("{project_dir}/{sample}/mapped/{sample}.sorted.bam", project_dir=config["project_dir"], sample=samples.index)
    output:
        expand("{project_dir}/intervals/SUMMARY/SELECTED_GENES_UNSORTED.bed", project_dir=config["project_dir"], sample=samples.index),
        expand("{project_dir}/intervals/{gene}.bed", project_dir=config["project_dir"], gene=selected_genes)
    conda:
        "../envs/intervals.yaml"
    script:
        "../scripts/intervals.py"

rule sort_selected_intervals:
    input:
        "{project_dir}/intervals/SUMMARY/SELECTED_GENES_UNSORTED.bed"
    output:
        "{project_dir}/intervals/SUMMARY/SELECTED_GENES.bed"
    shell:
        "sort -V -k 1,1 -k 2,2n -o {output} {input} \n"
        "rm {input}"


if _platform == "darwin" or config["variant_tool"] == "gatk":
    
    rule recalibrate_base_qualities:
        input:
            bam=get_recal_input(),
            bai=get_recal_input(bai=True),
            ref=rules.get_genome.output["fasta"],
            fasta_dict=rules.get_genome.output["fasta_dict"],
            known=rules.get_dbsnp.output["dbsnp"],
            tbi=rules.get_dbsnp.output["dbsnp_tbi"],
            regions="{project_dir}/intervals/SUMMARY/SELECTED_GENES.bed"
        output:
            recal_table="{project_dir}/{sample}/recal/{sample}.grp"
        params:
            extra=get_recal_params
        log:
            "{project_dir}/{sample}/logs/gatk/recal/{sample}.log"
        wrapper:
            "0.64.0/bio/gatk/baserecalibrator"

    rule gatk_applybqsr:
        input:
            bam=get_recal_input(),
            bai=get_recal_input(bai=True),
            ref=rules.get_genome.output["fasta"],
            fasta_dict=rules.get_genome.output["fasta_dict"],
            recal_table="{project_dir}/{sample}/recal/{sample}.grp"
        output:
            bam="{project_dir}/{sample}/applybqsr/{sample}.bam"
        log:
            "{project_dir}/{sample}/logs/gatk/applybqsr/{sample}.log"
        params:
            extra="",  # optional
            java_opts="", # optional
        wrapper:
            "0.64.0/bio/gatk/applybqsr"


rule samtools_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    wrapper:
        "0.64.0/bio/samtools/index"

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
        if config["processing"].get("restrict-regions"):
            # case 3: need an index because random access is required
            f += ".bai"
            return f
        else:
            # case 4: no index needed
            return []
    else:
        return f

def get_regions_param(regions=config["processing"].get("restrict-regions"), default=""):
    if regions:
        params = "--intervals '{}' ".format(regions)
        padding = config["processing"].get("region-padding")
        if padding:
            params += "--interval-padding {}".format(padding)
        return params
    return default


def get_regions_param(regions=config["processing"].get("restrict-regions"), default=""):
    if regions:
        params = "--regions '{}' ".format(regions)
        return params
    return default

##### END OF GATK SPECIFIC FUNCTIONS #####


if _platform == "darwin":
    rule map_reads:
        input:
            reads=get_map_reads_input,
            idx=rules.get_genome.output
        output:
            "{project_dir}/{sample}/mapped/{sample}.sorted.bam"
        log:
            "{project_dir}/{sample}/logs/bwa_mem/{sample}.log"
        params:
            index=lambda w, input: input.idx[0],
            extra=get_read_group,
            sort="samtools",
            sort_order="coordinate"
        threads: 8
        wrapper:
            "0.64.0/bio/bwa/mem"

elif _platform == "linux" or _platform == "linux2":
    rule map_reads:
        input:
            reads=get_map_reads_input,
            idx=rules.get_genome.output
        output:
            "{project_dir}/{sample}/mapped/{sample}.sorted.bam"
        log:
            "{project_dir}/{sample}/logs/bwa_mem2/{sample}.log"
        params:
            index=lambda w, input: input.idx[0],
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
        config["params"]["picard"]["MarkDuplicates"] + "USE_JDK_DEFLATER=true USE_JDK_INFLATER=true"
    wrapper:
        "0.64.0/bio/picard/markduplicates"


if _platform == "darwin" or config["variant_tool"] == "gatk":
    rule recalibrate_base_qualities:
        input:
            bam=get_recal_input(),
            bai=get_recal_input(bai=True),
            ref=rules.get_genome.output[0],
            idx=rules.genome_dict.output[0],
            known=rules.remove_iupac_codes.output[0],
            tbi=rules.tabix_known_variants.output[0]
        output:
            recal_table="{project_dir}/{sample}/recal/{sample}.grp"
        params:
            extra=get_regions_param() + config["params"]["gatk"]["BaseRecalibrator"]
        log:
            "{project_dir}/{sample}/logs/gatk/recal/{sample}.log"
        wrapper:
            "0.64.0/bio/gatk/baserecalibrator"

    rule gatk_applybqsr:
        input:
            bam=get_recal_input(),
            bai=get_recal_input(bai=True),
            ref=rules.get_genome.output[0],
            dict=rules.genome_dict.output[0],
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

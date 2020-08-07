##### GATK SPECIFIC FUNCTIONS #####

def get_sample_bams(wildcards):
    """Get all aligned reads of given sample."""
    return expand("{project_dir}/{sample}/recal/{sample}.bam",
                  sample=wildcards.sample)

def get_regions_param(regions=config["processing"].get("restrict-regions"), default=""):
    if regions:
        params = "--intervals '{}' ".format(regions)
        padding = config["processing"].get("region-padding")
        if padding:
            params += "--interval-padding {}".format(padding)
        return params
    return default

def get_call_variants_params(wildcards, input):
    return (get_regions_param(regions=input.regions, default="--intervals {}".format(wildcards.contig)) +
            config["params"]["gatk"]["HaplotypeCaller"])

def get_contigs():
    with checkpoints.genome_faidx.get().output[0].open() as fai:
        return pd.read_table(fai, header=None, usecols=[0], squeeze=True, dtype=str)

##### END OF GATK SPECIFIC FUNCTIONS #####


##### DEEPVARIANT SPECIFIC FUNCTIONS #####

def get_deepvariant_input(bai=False):
    # case 1: no duplicate removal
    f = "{project_dir}/{sample}/mapped/{sample}.sorted.bam"
    if config["processing"]["remove-duplicates"]:
        # case 2: remove duplicates
        f = "{project_dir}/{sample}/dedup/{sample}.bam"

def get_regions_param(regions=config["processing"].get("restrict-regions"), default=""):
    if regions:
        params = "--regions '{}' ".format(regions)
        return params
    return default

##### END OF DEEPVARIANT SPECIFIC FUNCTIONS #####


if _platform == "darwin":

    if "restrict-regions" in config["processing"]:
        rule compose_regions:
            input:
                config["processing"]["restrict-regions"]
            output:
                "{project_dir}/{sample}/intervals/{contig}.regions.bed"
            conda:
                "../envs/bedops.yaml"
            shell:
                "bedextract {wildcards.contig} {input} > {output}"

    rule call_variants:
        input:
            bam=get_sample_bams(),
            ref=rules.get_genome.output[0],
            idx=rules.genome_dict.output[0],
            known=rules.remove_iupac_codes.output[0],
            tbi=rules.tabix_known_variants.output[0],
            regions=rules.compose_regions.output[0] if config["processing"].get("restrict-regions") else []
        output:
            gvcf=protected("{project_dir}/{sample}/called/{sample}.{contig}.g.vcf.gz")
        log:
            "{project_dir}/{sample}/logs/gatk/haplotypecaller/{sample}.{contig}.log"
        params:
            extra=get_call_variants_params()
        wrapper:
            "0.64.0/bio/gatk/haplotypecaller"


    rule combine_calls:
        input:
            ref=rules.get_genome.output[0],
            gvcfs=expand("{project_dir}/{sample}/called/{sample}.{{contig}}.g.vcf.gz", sample=samples.index)
        output:
            gvcf="{project_dir}/{sample}/called/all.{contig}.g.vcf.gz"
        log:
            "{project_dir}/{sample}/logs/gatk/combinegvcfs.{contig}.log"
        wrapper:
            "0.64.0/bio/gatk/combinegvcfs"


    rule genotype_variants:
        input:
            ref=rules.get_genome.output[0],
            gvcf="{project_dir}/{sample}/called/all.{contig}.g.vcf.gz"
        output:
            vcf=temp("{project_dir}/{sample}/genotyped/all.{contig}.vcf.gz")
        params:
            extra=config["params"]["gatk"]["GenotypeGVCFs"]
        log:
            "{project_dir}/{sample}/logs/gatk/genotypegvcfs.{contig}.log"
        wrapper:
            "0.64.0/bio/gatk/genotypegvcfs"


    rule merge_variants:
        input:
            vcfs=lambda w: expand("{project_dir}/{sample}/genotyped/all.{contig}.vcf.gz", contig=get_contigs()),
        output:
            vcf="{project_dir}/{sample}/genotyped/all.vcf.gz"
        log:
            "{project_dir}/{sample}/logs/picard/merge-genotyped.log"
        wrapper:
            "0.64.0/bio/picard/mergevcfs"


elif _platform == "linux" or _platform == "linux2":

    rule deepvariant:
        input:
            bam=get_deepvariant_input(),
            ref=rules.get_genome.output[0]
        output:
            vcf="{project_dir}/{sample}/final_vcf/{sample}.vcf.gz"
        params:
            model=config["params"]["deepvariant"]["model"],
            extra=get_regions_param()
        threads: 2
        log:
            "{project_dir}/{sample}/logs/deepvariant/{sample}/stdout.log"
        wrapper:
            "0.64.0/bio/deepvariant"
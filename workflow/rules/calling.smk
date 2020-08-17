##### GATK SPECIFIC FUNCTIONS #####

def get_gatk_regions_param(regions="{project_dir}/intervals/{gene}.bed", default=""):
    if regions:
        params = "--intervals '{}' ".format(regions)
        padding = config["processing"].get("region-padding")
        if padding:
            params += "--interval-padding {}".format(padding)
        return params
    return default

def get_call_variants_params(wildcards, input):
    return (get_gatk_regions_param(regions=input.regions, default="--intervals {}".format(wildcards.gene)) +
            config["params"]["gatk"]["HaplotypeCaller"])

##### END OF GATK SPECIFIC FUNCTIONS #####


##### DEEPVARIANT SPECIFIC FUNCTIONS #####

def get_deepvariant_regions_param(regions="{project_dir}/intervals/{gene}.bed", default=""):
    if regions:
        params = "--regions '{}' ".format(regions)
        return params
    return default

def get_deepvariant_params(wildcards, input):
    return (get_deepvariant_regions_param(regions=input.regions, default="--regions {}".format(wildcards.gene)) +
            config["params"]["deepvariant"]["extra"])

##### END OF DEEPVARIANT SPECIFIC FUNCTIONS #####

##### COMMON FUNCTIONS #####

def get_sample_bams(wildcards):
    """
    Get the appropriate aligned bams of a given sample.
    """
    if _platform == "darwin" or config["variant_tool"] == "gatk":
        f = "{project_dir}/{sample}/applybqsr/{sample}.bam"
        return f
    elif (_platform == "linux" or _platform == "linux2") and config["variant_tool"] == "deepvariant":
        # Case 1: no duplicate removal
        f = "{project_dir}/{sample}/mapped/{sample}.sorted.bam"
        if config["processing"]["remove-duplicates"]:
            # case 2: remove duplicates
            f = "{project_dir}/{sample}/dedup/{sample}.bam"
        return f

##### END OF COMMON FUNCTIONS #####


if _platform == "darwin" or config["variant_tool"] == "gatk":

    rule call_variants:
        input:
            bam=get_sample_bams,
            ref=rules.get_genome.output["fasta"],
            idx=rules.get_genome.output["fasta_dict"],
            known=rules.get_dbsnp.output["dbsnp"],
            tbi=rules.get_dbsnp.output["dbsnp_tbi"],
            regions="{project_dir}/intervals/{gene}.bed"
        output:
            gvcf="{project_dir}/{sample}/genes/{gene}/gvcf/{gene}.g.vcf.gz"
        log:
            "{project_dir}/{sample}/logs/gatk/haplotypecaller/{gene}/{sample}.{gene}.log"
        params:
            extra=get_call_variants_params
        wrapper:
            "0.64.0/bio/gatk/haplotypecaller"


    rule genotype_variants:
        input:
            ref=rules.get_genome.output["fasta"],
            gvcf="{project_dir}/{sample}/genes/{gene}/gvcf/{gene}.g.vcf.gz"
        output:
            vcf="{project_dir}/{sample}/genes/{gene}/vcf/{gene}.vcf.gz"
        params:
            extra=config["params"]["gatk"]["GenotypeGVCFs"]
        log:
            "{project_dir}/{sample}/logs/gatk/genotypegvcfs/{gene}/{sample}.{gene}.log"
        wrapper:
            "0.64.0/bio/gatk/genotypegvcfs"


elif (_platform == "linux" or _platform == "linux2") and config["variant_tool"] == "deepvariant":

    rule deepvariant:
        input:
            bam=get_sample_bams,
            ref=rules.get_genome.output["fasta"],
            regions="{project_dir}/intervals/{gene}.bed"
        output:
            vcf="{project_dir}/{sample}/genes/{gene}/final_vcf/{gene}.vcf.gz"
        params:
            model=config["params"]["deepvariant"]["model"],
            extra=get_deepvariant_params
        threads: 2
        log:
            "{project_dir}/{sample}/logs/deepvariant/{gene}/{sample}.{gene}.log"
        wrapper:
            "0.64.0/bio/deepvariant"
def get_sample_bams(wildcards):
    """Get all aligned reads of given sample."""
    return expand("recal/{sample}/{sample}.bam",
                  sample=wildcards.sample)


def get_call_variants_params(wildcards, input):
    return (get_regions_param(regions=input.regions, default="--intervals {}".format(wildcards.contig)) +
            config["params"]["gatk"]["HaplotypeCaller"])

def get_contigs():
    with checkpoints.genome_faidx.get().output[0].open() as fai:
        return pd.read_table(fai, header=None, usecols=[0], squeeze=True, dtype=str)


if "restrict-regions" in config["processing"]:
    rule compose_regions:
        input:
            config["processing"]["restrict-regions"]
        output:
            "intervals/{contig}.regions.bed"
        conda:
            "../envs/bedops.yaml"
        shell:
            "bedextract {wildcards.contig} {input} > {output}"


rule call_variants:
    input:
        bam=get_sample_bams,
        ref=rules.get_genome.output[0],
        idx=rules.genome_dict.output[0],
        known=rules.remove_iupac_codes.output[0],
        tbi=rules.tabix_known_variants.output[0],
        regions="intervals/{contig}.regions.bed" if config["processing"].get("restrict-regions") else []
    output:
        gvcf=protected("called/{sample}/{sample}.{contig}.g.vcf.gz")
    log:
        "logs/gatk/haplotypecaller/{sample}/{sample}.{contig}.log"
    params:
        extra=get_call_variants_params
    wrapper:
        "0.64.0/bio/gatk/haplotypecaller"


rule combine_calls:
    input:
        ref=rules.get_genome.output[0],
        gvcfs=expand("called/{sample}/{sample}.{{contig}}.g.vcf.gz", sample=samples.index)
    output:
        gvcf="called/{sample}/all.{contig}.g.vcf.gz"
    log:
        "logs/gatk/{sample}/combinegvcfs.{contig}.log"
    wrapper:
        "0.64.0/bio/gatk/combinegvcfs"


rule genotype_variants:
    input:
        ref=rules.get_genome.output[0],
        gvcf="called/{sample}/all.{contig}.g.vcf.gz"
    output:
        vcf=temp("genotyped/{sample}/all.{contig}.vcf.gz")
    params:
        extra=config["params"]["gatk"]["GenotypeGVCFs"]
    log:
        "logs/gatk/{sample}/genotypegvcfs.{contig}.log"
    wrapper:
        "0.64.0/bio/gatk/genotypegvcfs"


rule merge_variants:
    input:
        vcfs=lambda w: expand("genotyped/{sample}/all.{contig}.vcf.gz", contig=get_contigs()),
    output:
        vcf="genotyped/{sample}/all.vcf.gz"
    log:
        "logs/picard/{sample}/merge-genotyped.log"
    wrapper:
        "0.64.0/bio/picard/mergevcfs"
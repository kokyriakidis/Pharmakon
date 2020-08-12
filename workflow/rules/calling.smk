##### GATK SPECIFIC FUNCTIONS #####

def get_gatk_regions_param(regions=config["processing"].get("restrict-regions"), default=""):
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

def get_deepvariant_regions_param(regions=config["processing"].get("restrict-regions"), default=""):
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


def read_selected_genes_table(fn: str) -> List[str]:
    """
    Read selected genes file.

    Returns:
        list[str]: A list of selected gene names.

    Args:
        fn (str): Selected genes file.
    """
    selected_genes = []
    with open(fn) as f:
        for gene in f: 
            gene = gene.strip()
            selected_genes.append(gene)

    return selected_genes


project_dir=config["project_dir"]

samples = pd.read_table(config["samples"]).set_index("sample_name", drop=False)
selected_samples = []
for sample in samples["sample_name"]:
    sample = sample.strip()
    selected_samples.append(sample)


#for sample in samples["sample_name"]:
selected_genes_table = f"{project_dir}/intervals/SUMMARY/SELECTED_GENES.txt"
selected_genes = read_selected_genes_table(selected_genes_table)



if _platform == "darwin" or config["variant_tool"] == "gatk":

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
            bam=get_sample_bams,
            ref=rules.get_genome.output[0],
            idx=rules.genome_dict.output[0],
            known=rules.remove_iupac_codes.output[0],
            tbi=rules.tabix_known_variants.output[0],
            regions=rules.compose_regions.output[0] if config["processing"].get("restrict-regions") else "{project_dir}/intervals/{gene}.bed"
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
            ref=rules.get_genome.output[0],
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
            ref=rules.get_genome.output[0]
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
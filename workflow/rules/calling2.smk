##### GATK SPECIFIC FUNCTIONS #####

# contigs in reference genome
def get_contigs():
    with checkpoints.genome_faidx.get().output[0].open() as fai:
        return pd.read_table(fai, header=None, usecols=[0], squeeze=True, dtype=str)

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

##### END OF GATK SPECIFIC FUNCTIONS #####


##### DEEPVARIANT SPECIFIC FUNCTIONS #####

def get_regions_param(regions=config["processing"].get("restrict-regions"), default=""):
    if regions:
        params = "--regions '{}' ".format(regions)
        return params
    return default

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

#from os.path import join

#DIR = expand("{project_dir}/{sample}", project_dir=config["project_dir"], sample=samples.index)

#DIR = config["project_dir"]
#SAMPLES = expand("/{sample}", sample=samples.index)
#for SAMPLE in SAMPLES:
#    GENES, = glob_wildcards(join(DIR, SAMPLE, "/intervals/{gene}.bed"))



#rule:
#	input: expand('{project_dir}/{sample}/{gene}.bed', project_dir=DIR, sample=SAMPLE ,gene=GENE)

#rule create_output:
#	input: lambda wildcards: [os.path.join(DIR[i], x + '.clean.vcf') for i,x in enumerate(NAME) if x == wildcards.name]
#	output: '{prefix}/{name}.output.vcf'
#	shell: 'cp {input} {output}'
#rule get_gene_intervals:
#    output:
#        expand("{project_dir}/{sample}/intervals/{gene}.bed", project_dir=config["project_dir"], sample=samples.index, gene=gene)
#    conda:
#        "../envs/intervals.yaml"
#    shell:
#        "python ../scripts/gene_intervals.py"




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
selected_genes_table = f"{project_dir}/intervals/SUMMARY/SELLECTED_GENES.txt"
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


    #rule regions:
    #    output:
    #        expand("{project_dir}/{sample}/intervals/{gene}.bed", project_dir=config["project_dir"], sample=samples.index, gene=selected_genes)
    #    run:
    #        project_dir=config["project_dir"]
    #        samples = pd.read_table(config["samples"]).set_index("sample_name", drop=False)
    #        gene_table = get_gene_table()
    #        for sample in samples["sample_name"]: 
    #            for gene in selected_genes:
    #                target_region =  "{chromosome}\t{start}\t{end}".format(chromosome=gene_table[gene]["chr"].replace("chr", ""), start=gene_table[gene]["hg38_start"], end=gene_table[gene]["hg38_end"])
    #                with open(
    #                    f"{project_dir}/{sample}/intervals/{gene}.bed", "w"
    #                ) as f:
    #                    f.write(target_region)




        rule call_variants:
            input:
                bam=get_sample_bams,
                ref=rules.get_genome.output[0],
                idx=rules.genome_dict.output[0],
                known=rules.remove_iupac_codes.output[0],
                tbi=rules.tabix_known_variants.output[0],
                #regions=rules.regions.output[0] if config["processing"].get("restrict-regions") else []
                regions=rules.compose_regions.output[0] if config["processing"].get("restrict-regions") else "{project_dir}/{sample}/intervals/{selected_genes}.bed"
            output:
                #gvcf="{project_dir}/{sample}/genes/{gene}/{sample}.{gene}.g.vcf.gz"
                gvcf=expand("{project_dir}/{sample}/genes/{selected_genes}/{sample}.{selected_genes}.g.vcf.gz", project_dir=config["project_dir"], sample=selected_samples, 
)

            log:
                #"{project_dir}/{sample}/logs/gatk/haplotypecaller/{sample}.{contig}.log"
            params:
                #extra=get_call_variants_params
            wrapper:
                "0.64.0/bio/gatk/haplotypecaller"


        rule combine_calls:
            input:
                ref=rules.get_genome.output[0],
                gvcfs=expand("{project_dir}/{sample}/called/{sample}.{{contig}}.g.vcf.gz", sample=samples.index, project_dir=config["project_dir"])
            output:
                gvcf="{project_dir}/{sample}/called/all.{contig}.g.vcf.gz"
            log:
                "{project_dir}/{sample}/logs/gatk/combinegvcfs/combinegvcfs.{contig}.log"
            wrapper:
                "0.64.0/bio/gatk/combinegvcfs"


        rule genotype_variants:
            input:
                ref=rules.get_genome.output[0],
                gvcf="{project_dir}/{sample}/called/all.{contig}.g.vcf.gz"
            output:
                vcf="{project_dir}/{sample}/genotyped/all.{contig}.vcf.gz"
            params:
                extra=config["params"]["gatk"]["GenotypeGVCFs"]
            log:
                "{project_dir}/{sample}/logs/gatk/genotypegvcfs/genotypegvcfs.{contig}.log"
            wrapper:
                "0.64.0/bio/gatk/genotypegvcfs"


        rule merge_variants:
            input:
                vcfs=lambda w: expand("{project_dir}/{sample}/genotyped/all.{contig}.vcf.gz", project_dir=config["project_dir"], sample=samples.index, contig=get_contigs()),
            output:
                vcf="{project_dir}/{sample}/genotyped/all.vcf.gz"
            log:
                "{project_dir}/{sample}/logs/picard/merge-genotyped.log"
            wrapper:
                "0.64.0/bio/picard/mergevcfs"


    elif (_platform == "linux" or _platform == "linux2") and config["variant_tool"] == "deepvariant":

        rule deepvariant:
            input:
                bam=get_sample_bams,
                ref=rules.get_genome.output[0]
            output:
                vcf="{project_dir}/{sample}/final_vcf/{sample}.vcf.gz"
            params:
                model=config["params"]["deepvariant"]["model"],
                extra=get_regions_param
            threads: 2
            log:
                "{project_dir}/{sample}/logs/deepvariant/{sample}/stdout.log"
            wrapper:
                "0.64.0/bio/deepvariant"
####################################

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


vartype=["snvs", "indels"]
filtertype="recalibrated" if config["filtering"]["vqsr"] else "hardfiltered"
######################################



if _platform == "darwin" or config["variant_tool"] == "gatk":

    def get_vartype_arg(wildcards):
        return "--select-type-to-include {}".format(
            "SNP" if wildcards.vartype == "snvs" else "INDEL")


    rule select_calls:
        input:
            ref=rules.get_genome.output["fasta"],
            vcf="{project_dir}/{sample}/genes/{gene}/vcf/{gene}.vcf.gz"
        output:
            vcf="{project_dir}/{sample}/genes/{gene}/filtered_vcf/{gene}.{vartype}.vcf.gz"
        params:
            extra=get_vartype_arg
        log:
            "{project_dir}/{sample}/logs/gatk/selectvariants/{gene}/{sample}.{gene}.{vartype}.log"
        wrapper:
            "0.59.0/bio/gatk/selectvariants"


    def get_filter(wildcards):
        return {
            "snv-hard-filter":
            config["filtering"]["hard"][wildcards.vartype]}


    rule hard_filter_calls:
        input:
            ref=rules.get_genome.output["fasta"],
            vcf="{project_dir}/{sample}/genes/{gene}/filtered_vcf/{gene}.{vartype}.vcf.gz"
        output:
            vcf="{project_dir}/{sample}/genes/{gene}/filtered_vcf/{gene}.{vartype}.hardfiltered.vcf.gz"
        params:
            filters=get_filter
        log:
           "{project_dir}/{sample}/logs/gatk/variantfiltration/{gene}/{sample}.{gene}.{vartype}.log"
        wrapper:
            "0.59.2/bio/gatk/variantfiltration"


    rule recalibrate_calls:
        input:
            vcf="{project_dir}/{sample}/genes/{gene}/filtered_vcf/{gene}.{vartype}.vcf.gz"
        output:
            vcf="{project_dir}/{sample}/genes/{gene}/filtered_vcf/{gene}.{vartype}.recalibrated.vcf.gz"
        params:
            extra=config["params"]["gatk"]["VariantRecalibrator"]
        log:
            "{project_dir}/{sample}/logs/gatk/variantrecalibrator/{gene}/{sample}.{gene}.{vartype}.log"
        wrapper:
            "0.59.2/bio/gatk/variantrecalibrator"


    rule merge_calls:
        input:
            expand("{project_dir}/{sample}/genes/{gene}/filtered_vcf/{gene}.{vartype}.{filtertype}.vcf.gz", project_dir=config["project_dir"], sample=samples.index, gene=selected_genes, vartype=["snvs", "indels"], filtertype="recalibrated" if config["filtering"]["vqsr"] else "hardfiltered")
        output:
            expand("{project_dir}/{sample}/genes/{gene}/final_vcf/{gene}.vcf.gz", project_dir=config["project_dir"], sample=samples.index, gene=selected_genes)
        conda:
            "../envs/mergevcfs.yaml"
        script:
            "../scripts/mergevcfs.py"







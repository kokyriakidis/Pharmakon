def get_gdf_input(wildcards):
    if _platform == "darwin" or config["variant_tool"] == "gatk":
        f = "{project_dir}/{sample}/applybqsr/{sample}.bam"
    else:
        # case 1: no duplicate removal
        f = "{project_dir}/{sample}/mapped/{sample}.sorted.bam"
        if config["processing"]["remove-duplicates"]:
            # case 2: remove duplicates
            f = "{project_dir}/{sample}/dedup/{sample}.bam"
    return f


rule get_gdf:
    input:
        expand("{project_dir}/{sample}/genes/{gene}/final_vcf/{gene}.vcf.gz", project_dir=config["project_dir"], sample=samples.index, gene=selected_genes)
    output:
        expand("{project_dir}/{sample}/genes/{gene}/gdf/{gene}.gdf", project_dir=config["project_dir"], sample=samples.index, gene=selected_genes)
    conda:
        "../envs/gdf.yaml"
    script:
        "../scripts/bam2gdf.py"

rule get_stargazer:
    output:
        "stargazer/setup.py"
    conda:
        "../envs/stargazer.yaml"
    shell:
        """
        git clone -b "v1.2.1" --single-branch https://github.com/sbslee/stargazer.git
        cd stargazer
        python setup.py install
        """

rule run_stargazer:
    input:
        "stargazer/setup.py",
        expand("{project_dir}/{sample}/genes/{gene}/gdf/{gene}.gdf", project_dir=config["project_dir"], sample=samples.index, gene=selected_genes)
    output:
        expand("{project_dir}/{sample}/genes/{gene}/stargazer/genotype.txt", project_dir=config["project_dir"], sample=samples.index, gene=selected_genes)
    conda:
        "../envs/stargazer.yaml"
    script:
        "../scripts/stargazer.py"


rule merge_stargazer_genotypes:
    input:
        expand("{project_dir}/{sample}/genes/{gene}/stargazer/genotype.txt", project_dir=config["project_dir"], sample=samples.index, gene=selected_genes)
    output:
        expand("{project_dir}/{sample}/genotypes/genotypes.txt", project_dir=config["project_dir"], sample=samples.index)
    script:
        "../scripts/merge_genotypes.py"



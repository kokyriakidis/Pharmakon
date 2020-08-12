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

#rule get_gdf:
#    output:
#        expand("{project_dir}/{sample}/genes/{gene}/gdf/{gene}.gdf", project_dir=config["project_dir"], sample=samples.index, gene=selected_genes)
#	conda:
#		"../envs/gdf.yaml"
#    script:
#        "../scripts/bam2gdf2.py"
#    #log:
#    #    "{project_dir}/{sample}/logs/gdf/{gene}.gdf.log"

rule get_gdf:
    output:
        expand("{project_dir}/{sample}/genes/{gene}/gdf/{gene}.gdf", project_dir=config["project_dir"], sample=samples.index, gene=selected_genes)
    conda:
        "../envs/gdf.yaml"
    log:
        expand("{project_dir}/{sample}/logs/gdf/{gene}.gdf.log", project_dir=config["project_dir"], sample=samples.index, gene=selected_genes)
    script:
        "../scripts/bam2gdf2.py"

#rule stargazer:
#    input:
#        stargazer_tool=config["stargazer"]["path"] + "/stargazer"
#        data_type="wgs"
#        genome_build="hg38" if config["ref"]["build"] == "GRCh38" else "hg19",
#        target_gene=config["params"]["stargazer"]["target_gene"],
#        vcf_file="{project_dir}/{sample}/final_vcf/{sample}.vcf.gz",
#        control_gene=config["params"]["stargazer"]["control_gene"],
#        gdf_file="{project_dir}/{sample}/gdf/{sample}.gdf"
#        project_dir="{project_dir}/{sample}/stargazer"
#    output:
#        output_file="{project_dir}/{sample}/stargazer/{sample}_genotype.txt"
#    threads: 2
#    log:
#        "{project_dir}/{sample}/logs/stargazer/stargazer.log"
#    shell: 
#        "{input.stargazer_tool} {input.data_type} {input.genome_build} {input.target_gene} {input.vcf_file} {input.project_dir} --cg {input.control_gene} --gdf {input.gdf_file} --plot"

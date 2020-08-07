rule get_gdf:
    input:
        bam_file="dedup/{sample}/{sample}.bam",
        genome_build=config["ref"]["build"],
        target_gene=config["params"]["stargazer"]["target_gene"],
        control_gene=config["params"]["stargazer"]["control_gene"],
    output:
        output_file="stargazer/{sample}/gdf/{sample}.gdf"
    threads: 2
    log:
        "logs/stargazer/{sample}/get_gdf.log"
    script:
        "scripts/bam2gdf.py"


rule stargazer:
    input:
        stargazer_tool=config["stargazer"]["path"] + "/stargazer"
        data_type="wgs"
        genome_build="hg38" if config["ref"]["build"] == "GRCh38" else "hg19",
        target_gene=config["params"]["stargazer"]["target_gene"],
        vcf_file="calls/{sample}/{sample}.vcf.gz",
        control_gene=config["params"]["stargazer"]["control_gene"],
        gdf_file="stargazer/{sample}/gdf/{sample}.gdf"
        project_dir="stargazer/{sample}/stargazer"
    output:
        output_file="stargazer/{sample}/stargazer/{sample}_genotype.txt"
    threads: 2
    log:
        "logs/stargazer/{sample}/stargazer.log"
    shell: 
        "{input.stargazer_tool} {input.data_type} {input.genome_build} {input.target_gene} {input.vcf_file} {input.project_dir} --cg {input.control_gene} --gdf {input.gdf_file} --plot"

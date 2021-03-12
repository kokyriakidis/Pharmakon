rule gdf__calculate_read_depth:
    """"
    Create a GDF (GATK DepthOfCoverage Format) file for Stargazer from BAM files
    by computing read depth.
    :input bam: Read BAM file
    :input target_gene: Name of the target gene
    :input control_gene: Name of a preselected control gene
    :input genome_build: Build of the reference genome assembly
    :output gdf: Path to the output GDF file
    """"
    input:
        bam          = f"{OUTDIR}/{{sample}}/{{mapper}}/{{sample}}.sorted.bam",
        target_gene  = f"{target_gene}",
        control_gene = f"{control_gene}",
        genome_build = config["build"]
    output: 
        gdf = f"{OUTDIR}/{{sample}}/stargazer/gdf/{{target_gene}}/{{sample}}.{{target_gene}}.gdf"
    log: 
        f"{OUTDIR}/{{sample}}/logs/stargazer/gdf/{{target_gene}}/{{sample}}.{{target_gene}}.gdf.log"
    benchmark: 
        f"{OUTDIR}/{{sample}}/logs/stargazer/gdf/{{target_gene}}/{{sample}}.{{target_gene}}.gdf.benchmark"
    conda: 
        "envs/pypgx.yaml"
    shell: 
        """
        (pypgx calculate-read-depth -t {input.target_gene} -c {input.control_gene} -i {input.bam} -o {output.gdf}) > {log} 2>&1
        """


rule stargazer__call_star_alleles:
    """"
    Stargazer is a bioinformatics tool for calling star alleles (haplotypes) in PGx genes
    :input stargazer_path: Path to Stargazer folder
    :input vcf: VCF file with called germline variants
    :input gdf: GDF file of a PGx gene
    :input target_gene: Name of the target gene
    :input control_gene: Name of a preselected control gene
    :input datatype: Targeted sequencing (ts) or Whole Genome Sequencing (wgs)
    :output output_dir: Output files directory
    :output output_prefix: Output filename prefix
    :output genotype: Star alleles (haplotypes) in a text format
    """"
    input:
        stargazer_path = config["stargazer_path"],
        vcf            = f"{OUTDIR}/{{sample}}/{{caller}}/{{sample}}.{{caller}}.vcf.gz",
        gdf            = f"{OUTDIR}/{{sample}}/stargazer/gdf/{{target_gene}}/{{sample}}.{{target_gene}}.gdf",
        target_gene    = f"{target_gene}",
        control_gene   = f"{control_gene}",
        datatype       = "ts" if config["datatype"] == "ILLUMINA-WES" else "wgs"
    output:
        output_dir    = f"{OUTDIR}/{{sample}}/stargazer/genotype/{{target_gene}}",
        output_prefix = f"{sample}.{target_gene}.genotype",
        genotype      = f"{OUTDIR}/{{sample}}/stargazer/genotype/{{target_gene}}/{{sample}}.{{target_gene}}.genotype.txt"
    log: 
        f"{OUTDIR}/{{sample}}/logs/stargazer/genotype/{{target_gene}}/{{sample}}.{{target_gene}}.genotype.log"
    benchmark: 
        f"{OUTDIR}/{{sample}}/logs/stargazer/genotype/{{target_gene}}/{{sample}}.{{target_gene}}.genotype.benchmark"
    conda:
        "../envs/stargazer.yaml"
    shell: 
        """
        (python {input.stargazer_path}/stargazer.py genotype \
        --output_dir {output.output_dir} \
        --output_prefix {output.output_prefix} \
        -d {input.datatype} \
        -t {input.target_gene} \
        --vcf {input.vcf} \
        -c {input.control_gene} \
        --detail \
        --gdf GDF) > {log} 2>&1
        """


rule combine__merge_stargazer_genotypes:
    input:
        expand(f"{OUTDIR}/{{sample}}/stargazer/genotype/{{target_gene}}/{{sample}}.{{target_gene}}.genotype.txt", sample=wildcards.sample, target_gene=TARGET_GENES)

    output:
        expand(f"{OUTDIR}/{{sample}}/stargazer/{{sample}}.merged.genotypes.txt", sample=wildcards.sample)
    script:
        "../scripts/merge_genotypes.py"
    run:
        for sample in config["samples"]:
            random_selected_gene = TARGET_GENES[0]
            shell("head -n1 f'{OUTDIR}/{{sample}}/stargazer/genotype/{random_selected_gene}/stargazer/genotype.txt' > {project_dir}/{sample}/genotypes/genotypes.txt\n")
            for gene in selected_genes:
                shell(
                    "tail -n+2 {project_dir}/{sample}/genes/{gene}/stargazer/genotype.txt >> {project_dir}/{sample}/genotypes/genotypes.txt\n"
                    )



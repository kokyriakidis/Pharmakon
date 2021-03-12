from multiprocessing import cpu_count

#def get_bams(wildcards):
#    """Get all aligned reads of given sample."""
#    if config["datatype"] == "ILLUMINA-WGS" or config["datatype"] == "ILLUMINA-WES" or config["datatype"] == "NANOPORE":
#        return expand(f"{OUTDIR}/{{sample}}/minimap2/{{sample}}.sorted.bam", sample=wildcards.sample)
#    elif config["datatype"] == "PACBIO":
#        return expand(f"{OUTDIR}/{{sample}}/pbmm2/{{sample}}.sorted.bam", sample=wildcards.sample)


rule fastqc__quality_control:
    """
    Aggregate quality control results from all FastQC reports and generate summary HTML table
    :input reads: Sequenced reads in fastq format
    :output html: Quality report in HTML format
    :output zip: Aggregate quality control results
    """
    input:
        reads = get_samples_data
    output:
        html = f"{OUTDIR}/{{sample}}/qc/fastqc/{{sample}}_fastqc.html",
        zip  = f"{OUTDIR}/{{sample}}/qc/fastqc/{{sample}}_fastqc.zip"
    log:
        f"{OUTDIR}/{{sample}}/logs/fastqc/{{sample}}.fastqc.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/logs/fastqc/{{sample}}.fastqc.benchmark"
    threads:
        lambda cores: cpu_count() - 2
    wrapper:
        "0.72.0/bio/fastqc"


rule samtools__bam_stats:
    """
    Generate stats using samtools
    :input bams: Read alignment files in a binary format
    :output stats: Statistics from BAM files outputted in a text format
    """
    input:
        bams = f"{OUTDIR}/{{sample}}/{{mapper}}/{{sample}}.sorted.bam"
    output:
        stats = f"{OUTDIR}/{{sample}}/qc/samtools/{{sample}}.bam.stats.txt"
    log:
        f"{OUTDIR}/{{sample}}/logs/samtools/{{sample}}.bam.stats.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/logs/samtools/{{sample}}.bam.stats.benchmark"
    wrapper:
        "0.72.0/bio/samtools/stats"


rule bcftools_stats:
    """
    Calculating VCF statistics using bcftools
    :input bams: Read alignment files in a binary format
    :output stats: Statistics from VCF files outputted in a text format
    """
    input: 
        vcf = f"{OUTDIR}/{{sample}}/{{caller}}/{{sample}}.{{caller}}.vcf.gz"
    output: 
        stats = f"{OUTDIR}/{{sample}}/qc/bcftools/{{sample}}.{{caller}}.vcf.stats.txt"
    log: 
        f"{OUTDIR}/{{sample}}/logs/bcftools/{{sample}}.{{caller}}.vcf.stats.log"
    benchmark: 
        f"{OUTDIR}/{{sample}}/logs/bcftools/{{sample}}.{{caller}}.vcf.stats.benchmark"
    params: 
        fasta   = f"--fasta-ref {config['fasta']}", 
        filters = "--apply-filters PASS", 
        sample  = f"-s {{sample}}"
    threads:
        lambda cores: cpu_count() - 2
    conda: 
        "envs/bcftools.yaml"
    message: "Executing {rule}: Calculating VCF statistics for {input}."
    shell: 
        """
        (bcftools stats --threads {threads} {params.fasta} {params.filters} {params.sample} {input.vcf} > {output.stats}) > {log} 2>&1
        """


rule multiqc:
    """
    Generate QC report using multiqc
    """
    input:
        expand([f"{OUTDIR}/{{sample}}/qc/samtools/{{sample}}.bam.stats.txt",
                f"{OUTDIR}/{{sample}}/qc/bcftools/{{sample}}.{{caller}}.vcf.stats.txt",
                f"{OUTDIR}/{{sample}}/qc/fastqc/{{sample}}_fastqc.zip"],
                sample=f"{{sample}}")
    output:
        expand(f"{OUTDIR}/{{sample}}/qc/multiqc/multiqc.html", sample=wildcards.sample)
    log:
        expand(f"{OUTDIR}/{{sample}}/logs/multiqc/multiqc.log", sample=wildcards.sample)
    benchmark:
        expand(f"{OUTDIR}/{{sample}}/logs/multiqc/multiqc.benchmark", sample=wildcards.sample)
    wrapper:
        "0.72.0/bio/multiqc"

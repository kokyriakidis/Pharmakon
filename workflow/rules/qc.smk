from multiprocessing import cpu_count

## ILLUMINA QC ##

rule fastqc__reads_quality_control:
    """
    Quality control results from Illumina reads and generation of summary HTML table
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
        f"{OUTDIR}/{{sample}}/qc/fastqc/logs/{{sample}}.fastqc.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/qc/fastqc/logs/{{sample}}.fastqc.benchmark"
    threads:
        lambda cores: cpu_count() - 2
    wrapper:
        "0.72.0/bio/fastqc"



## NANOPORE AND PACBIO QC ##

rule nanoplot__bam_statistics:
    """
    Calculate various statistics from a long read sequencing dataset in BAM format.
    Creates a statistical summary using NanoStats.
    Creates various plots with NanoPlot.
    Creates a html report based on the previously created plots.
    """
    input:
        bam = f"{OUTDIR}/{{sample}}/{MAPPER}/{{sample}}.{BUILD}.{PROVIDER}.sorted.bam",
        bai = f"{OUTDIR}/{{sample}}/{MAPPER}/{{sample}}.{BUILD}.{PROVIDER}.sorted.bam.bai"
    output:
        f"{OUTDIR}/{{sample}}/qc/nanoplot/{{sample}}.{BUILD}.{PROVIDER}.bam.NanoStats.txt",
        f"{OUTDIR}/{{sample}}/qc/nanoplot/{{sample}}.{BUILD}.{PROVIDER}.bam.NanoPlot-report.html"
    params:
        output_dir    = f"{OUTDIR}/{{sample}}/qc/nanoplot",
        output_prefix = f"{{sample}}.{BUILD}.{PROVIDER}.bam."
    benchmark:
        f"{OUTDIR}/{{sample}}/qc/nanoplot/logs/{{sample}}.{BUILD}.{PROVIDER}.bam.nanoplot.benchmark"
    threads:
        lambda cores: cpu_count() - 2
    conda: 
        "../envs/nanoplot.yaml"
    shell:
        """
        NanoPlot --bam {input.bam} --outdir {params.output_dir} --prefix {params.output_prefix} --threads {threads}
        mv {params.output_dir}/*.log {params.output_dir}/logs/
        """


## COMMON QC ##

rule samtools__bam_stats:
    """
    Generate stats using samtools
    :input bams: Read alignment files in a binary format
    :output stats: Statistics from BAM files outputted in a text format
    """
    input:
        bam = f"{OUTDIR}/{{sample}}/{MAPPER}/{{sample}}.{BUILD}.{PROVIDER}.sorted.bam",
        bai = f"{OUTDIR}/{{sample}}/{MAPPER}/{{sample}}.{BUILD}.{PROVIDER}.sorted.bam.bai"
    output:
        stats     = f"{OUTDIR}/{{sample}}/qc/samtools/{{sample}}.{BUILD}.{PROVIDER}.{MAPPER}.bam.stats.txt",
        idxstats  = f"{OUTDIR}/{{sample}}/qc/samtools/{{sample}}.{BUILD}.{PROVIDER}.{MAPPER}.bam.idxstats.txt",
        flagstats = f"{OUTDIR}/{{sample}}/qc/samtools/{{sample}}.{BUILD}.{PROVIDER}.{MAPPER}.bam.flagstats.txt"
    log:
        f"{OUTDIR}/{{sample}}/qc/samtools/logs/{{sample}}.{BUILD}.{PROVIDER}.{MAPPER}.bam.stats.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/qc/samtools/logs/{{sample}}.{BUILD}.{PROVIDER}.{MAPPER}.bam.stats.benchmark"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools idxstats {input.bam} > {output.idxstats}
        samtools stats {input.bam} > {output.stats}
        samtools flagstat {input.bam} > {output.flagstats}
        """

rule mosdepth__calculate_coverage:
    """
    Calculating aligned coverage depth using mosdepth.
    """
    input:
        bam = f"{OUTDIR}/{{sample}}/{MAPPER}/{{sample}}.{BUILD}.{PROVIDER}.sorted.bam",
        bai = f"{OUTDIR}/{{sample}}/{MAPPER}/{{sample}}.{BUILD}.{PROVIDER}.sorted.bam.bai"
    output:
        f"{OUTDIR}/{{sample}}/qc/mosdepth/{{sample}}.{BUILD}.{PROVIDER}.mosdepth.global.dist.txt",
        f"{OUTDIR}/{{sample}}/qc/mosdepth/{{sample}}.{BUILD}.{PROVIDER}.mosdepth.region.dist.txt",
        f"{OUTDIR}/{{sample}}/qc/mosdepth/{{sample}}.{BUILD}.{PROVIDER}.mosdepth.summary.txt",
        f"{OUTDIR}/{{sample}}/qc/mosdepth/{{sample}}.{BUILD}.{PROVIDER}.regions.bed.gz"
    log: 
        f"{OUTDIR}/{{sample}}/qc/mosdepth/logs/{{sample}}.{BUILD}.{PROVIDER}.mosdepth.log"
    benchmark: 
        f"{OUTDIR}/{{sample}}/qc/mosdepth/logs/{{sample}}.{BUILD}.{PROVIDER}.mosdepth.benchmark"
    params:
        by     = "500",
        prefix = f"{OUTDIR}/{{sample}}/qc/mosdepth/{{sample}}.{BUILD}.{PROVIDER}",
        extra  = "--no-per-base --use-median"
    threads: 4
    conda: 
        "../envs/mosdepth.yaml"
    shell:
        """
        (mosdepth \
            --threads {threads} --by {params.by} \
            {params.extra} {params.prefix} {input.bam}) > {log} 2>&1
        """

rule sex__infer_sex_from_coverage:
    """
    Inference of chromosomal sex from mosdepth coverage summary to check for sample swaps.
    """
    input: 
        f"{OUTDIR}/{{sample}}/qc/mosdepth/{{sample}}.{BUILD}.{PROVIDER}.mosdepth.summary.txt"
    output: 
        f"{OUTDIR}/{{sample}}/qc/mosdepth/{{sample}}.{BUILD}.{PROVIDER}.mosdepth.inferred_sex.txt"
    log: 
        f"{OUTDIR}/{{sample}}/qc/mosdepth/logs/{{sample}}.{BUILD}.{PROVIDER}.infer_sex_from_coverage.log"
    benchmark: 
        f"{OUTDIR}/{{sample}}/qc/mosdepth/logs/{{sample}}.{BUILD}.{PROVIDER}.infer_sex_from_coverage.benchmark"
    conda: 
        "../envs/pandas.yaml"
    shell: 
        """
        (python3 workflow/scripts/infer_sex_from_coverage.py {input} > {output}) > {log} 2>&1
        """


rule consistency__calculate_m2_ratio:
    """
    Calculation of chrM:chr2 ratio from mosdepth coverage summary to check for consistency between runs.
    """
    input: 
        f"{OUTDIR}/{{sample}}/qc/mosdepth/{{sample}}.{BUILD}.{PROVIDER}.mosdepth.summary.txt"
    output: 
        f"{OUTDIR}/{{sample}}/qc/mosdepth/{{sample}}.{BUILD}.{PROVIDER}.mosdepth.M2_ratio.txt"
    log: 
        f"{OUTDIR}/{{sample}}/qc/mosdepth/logs/{{sample}}.{BUILD}.{PROVIDER}.calculate_M2_ratio.log"
    benchmark: 
        f"{OUTDIR}/{{sample}}/qc/mosdepth/logs/{{sample}}.{BUILD}.{PROVIDER}.calculate_M2_ratio.benchmark"
    conda: 
        "../envs/pandas.yaml"
    shell: 
        """
        (python3 workflow/scripts/calculate_M2_ratio.py {input} > {output}) > {log} 2>&1
        """

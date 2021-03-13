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
        bam = f"{OUTDIR}/{{sample}}/{{mapper}}/{{sample}}.{{provider}}.sorted.bam",
        bai = f"{OUTDIR}/{{sample}}/{{mapper}}/{{sample}}.{{provider}}.sorted.bam.bai"
    output:
        stats     = f"{OUTDIR}/{{sample}}/qc/samtools/{{sample}}.{{mapper}}.{{provider}}.bam.stats.txt",
        idxstats  = f"{OUTDIR}/{{sample}}/qc/samtools/{{sample}}.{{mapper}}.{{provider}}.bam.idxstats.txt",
        flagstats = f"{OUTDIR}/{{sample}}/qc/samtools/{{sample}}.{{mapper}}.{{provider}}.bam.flagstats.txt"
    log:
        f"{OUTDIR}/{{sample}}/logs/samtools/{{sample}}.{{mapper}}.{{provider}}.bam.stats.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/logs/samtools/{{sample}}.{{mapper}}.{{provider}}.bam.stats.benchmark"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools idxstats {input.bam} > {output.idxstats}
        samtools stats {input.bam} > {output.stats}
        samtools flagstat {input.bam} > {output.flagstats}
        """

rule bcftools_stats:
    """
    Calculating VCF statistics using bcftools
    :input bams: Read alignment files in a binary format
    :output stats: Statistics from VCF files outputted in a text format
    """
    input: 
        vcf = f"{OUTDIR}/{{sample}}/{{caller}}/{{sample}}.{{provider}}.{{caller}}.vcf.gz"
    output: 
        stats = f"{OUTDIR}/{{sample}}/qc/bcftools/{{sample}}.{{provider}}.{{caller}}.vcf.stats.txt"
    log: 
        f"{OUTDIR}/{{sample}}/logs/bcftools/{{sample}}.{{provider}}.{{caller}}.vcf.stats.log"
    benchmark: 
        f"{OUTDIR}/{{sample}}/logs/bcftools/{{sample}}.{{provider}}.{{caller}}.vcf.stats.benchmark"
    params: 
        fasta   = f"--fasta-ref {config['fasta']}", 
        filters = "--apply-filters PASS", 
        sample  = f"-s {{sample}}"
    threads:
        lambda cores: cpu_count() - 2
    conda: 
        "../envs/bcftools.yaml"
    message: "Executing {rule}: Calculating VCF statistics for {input}."
    shell: 
        """
        (bcftools stats --threads {threads} {params.fasta} {params.filters} {params.sample} {input.vcf} > {output.stats}) > {log} 2>&1
        """

## PACBIO QC ##

rule smrtcell__fastq_stats:
    """
    Generate read length and quality stats
    """
    input: 
        get_samples_data
    output: 
        f"{OUTDIR}/{{sample}}/qc/smrtcell_stats/{{sample}}.read_length_and_quality.tsv"
    log: 
        f"{OUTDIR}/{{sample}}/logs/smrtcell_stats/{{sample}}.read_length_and_quality.log"
    benchmark: 
        f"{OUTDIR}/{{sample}}/logs/smrtcell_stats/{{sample}}.read_length_and_quality.benchmark"
    conda: 
        "../envs/smrtcell_stats.yaml"
    shell: 
        """
        (python3 workflow/scripts/extract_read_length_and_qual.py {input} > {output}) > {log} 2>&1
        """


rule smrtcell_summary_stats:
    """
    Summarize read length and quality stats
    """
    input: 
        f"{OUTDIR}/{{sample}}/qc/smrtcell_stats/{{sample}}.read_length_and_quality.tsv"
    output:
        rlsummary = f"{OUTDIR}/{{sample}}/qc/smrtcell_stats/{{sample}}.read_length_summary.tsv",
        rqsummary = f"{OUTDIR}/{{sample}}/qc/smrtcell_stats/{{sample}}.read_quality_summary.tsv"
    log: 
        f"{OUTDIR}/{{sample}}/logs/smrtcell_stats/{{sample}}.summary.log"
    benchmark: 
        f"{OUTDIR}/{{sample}}/logs/smrtcell_stats/{{sample}}.summary.tsv"
    conda: 
        "../envs/smrtcell_stats.yaml"
    shell:
        """
        (awk '{{ b=int($2/1000); b=(b>39?39:b); print 1000*b "\t" $2; }}' {input} |
            sort -k1,1g | datamash -g 1 count 1 sum 2 |
            awk 'BEGIN {{ for(i=0;i<=39;i++) {{ print 1000*i"\t0\t0"; }} }} {{ print; }}' |
            sort -k1,1g | datamash -g 1 sum 2 sum 3 > {output.rlsummary}) 2> {log}
        (awk '{{ print ($3>50?50:$3) "\t" $2; }}' {input} |
            sort -k1,1g | datamash -g 1 count 1 sum 2 |
            awk 'BEGIN {{ for(i=0;i<=60;i++) {{ print i"\t0\t0"; }} }} {{ print; }}' |
            sort -k1,1g | datamash -g 1 sum 2 sum 3 > {output.rqsummary}) 2>> {log}
        """

rule mosdepth__calculate_coverage:
    """
    Calculating coverage using mosdepth
    """
    input:
        bam = f"{OUTDIR}/{{sample}}/pbmm2/{{sample}}.pacbio.sorted.bam",
        bai = f"{OUTDIR}/{{sample}}/pbmm2/{{sample}}.pacbio.sorted.bam.bai"
    output:
        f"{OUTDIR}/{{sample}}/qc/mosdepth/{{sample}}.{config['build']}.mosdepth.global.dist.txt",
        f"{OUTDIR}/{{sample}}/qc/mosdepth/{{sample}}.{config['build']}.mosdepth.region.dist.txt",
        f"{OUTDIR}/{{sample}}/qc/mosdepth/{{sample}}.{config['build']}.mosdepth.summary.txt",
        f"{OUTDIR}/{{sample}}/qc/mosdepth/{{sample}}.{config['build']}.regions.bed.gz"
    log: 
        f"{OUTDIR}/{{sample}}/logs/mosdepth/{{sample}}.{config['build']}.mosdepth.log"
    benchmark: 
        f"{OUTDIR}/{{sample}}/logs/mosdepth/{{sample}}.{config['build']}.mosdepth.benchmark"
    params:
        by = "500",
        prefix = f"{OUTDIR}/{{sample}}/qc/mosdepth/{{sample}}.{config['build']}",
        extra = "--no-per-base --use-median"
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
    Inference of chromosomal sex from mosdepth coverage summary
    """
    input: 
        f"{OUTDIR}/{{sample}}/qc/mosdepth/{{sample}}.{config['build']}.mosdepth.summary.txt"
    output: 
        f"{OUTDIR}/{{sample}}/qc/mosdepth/{{sample}}.{config['build']}.mosdepth.inferred_sex.txt"
    log: 
        f"{OUTDIR}/{{sample}}/logs/quality_control/{{sample}}.{config['build']}.infer_sex_from_coverage.log"
    benchmark: 
        f"{OUTDIR}/{{sample}}/logs/quality_control/{{sample}}.{config['build']}.infer_sex_from_coverage.benchmark"
    conda: 
        "../envs/pandas.yaml"
    shell: 
        """
        (python3 workflow/scripts/infer_sex_from_coverage.py {input} > {output}) > {log} 2>&1
        """


rule consistency__calculate_m2_ratio:
    """
    Calculation of chrM:chr2 ratio from mosdepth coverage summary
    """
    input: 
        f"{OUTDIR}/{{sample}}/qc/mosdepth/{{sample}}.{config['build']}.mosdepth.summary.txt"
    output: 
        f"{OUTDIR}/{{sample}}/qc/mosdepth/{{sample}}.{config['build']}.mosdepth.M2_ratio.txt"
    log: 
        f"{OUTDIR}/{{sample}}/logs/quality_control/{{sample}}.{config['build']}.calculate_M2_ratio.log"
    benchmark: 
        f"{OUTDIR}/{{sample}}/logs/quality_control/{{sample}}.{config['build']}.calculate_M2_ratio.benchmark"
    conda: 
        "../envs/pandas.yaml"
    shell: 
        """
        (python3 workflow/scripts/calculate_M2_ratio.py {input} > {output}) > {log} 2>&1
        """


rule coverage__calculate_gc_coverage:
    """
    Calculation of GC coverage distribution from mosdepth coverage by region
    """
    input:
        mosdepth_regions = f"{OUTDIR}/{{sample}}/qc/mosdepth/{{sample}}.{config['build']}.regions.bed.gz",
        ref = config["fasta"]
    output: 
        f"{OUTDIR}/{{sample}}/qc/mosdepth/{{sample}}.{config['build']}.gc_coverage.summary.txt"
    log: 
        f"{OUTDIR}/{{sample}}/logs/quality_control/{{sample}}.{config['build']}.calculate_gc_coverage.log"
    benchmark: 
        f"{OUTDIR}/{{sample}}/logs/quality_control/{{sample}}.{config['build']}.calculate_gc_coverage.benchmark"
    conda: 
        "../envs/gc_coverage.yaml"
    shell:
        """
        (bedtools nuc -fi {input.ref} -bed {input.mosdepth_regions} \
            | awk '($11==0) {{ print 0.05*int($6/0.05) "\t" $4; }}' \
            | sort -k1,1g \
            | datamash -g1 q1 2 median 2 q3 2 count 2 \
            | tr '\t' ',' \
            | awk 'BEGIN {{ print "#gc_perc,q1,median,q3,count"; }} {{ print $0; }}' > {output}) > {log} 2>&1
        """

## NANOPORE ##

rule nanoplot__bam_statistics:
    """
    Calculate various statistics from a long read sequencing dataset in BAM format.
    Creates a statistical summary using NanoStats.
    Creates various plots with nanoplotter.
    Creates a html report based on the previously created plots.
    """
    input:
        bam           = f"{OUTDIR}/{{sample}}/minimap2/{{sample}}.nanopore.sorted.bam",
        bai           = f"{OUTDIR}/{{sample}}/minimap2/{{sample}}.nanopore.sorted.bam.bai"
    output:
        f"{OUTDIR}/{{sample}}/qc/nanoplot/{{sample}}.bam.NanoStats.txt",
        f"{OUTDIR}/{{sample}}/qc/nanoplot/{{sample}}.bam.NanoPlot-report.html"
    params:
        output_dir    = f"{OUTDIR}/{{sample}}/qc/nanoplot",
        output_prefix = f"{{sample}}.bam."
    log:
        f"{OUTDIR}/{{sample}}/logs/nanoplot/{{sample}}.bam.nanoplot.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/logs/nanoplot/{{sample}}.bam.nanoplot.benchmark"
    threads:
        lambda cores: cpu_count() - 2
    conda: 
        "../envs/nanoplot.yaml"
    shell:
        """
        (NanoPlot --bam {input.bam} --N50 --outdir {params.output_dir} --prefix {params.output_prefix} --threads {threads}) > {log} 2>&1
        """
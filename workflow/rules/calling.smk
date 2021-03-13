from multiprocessing import cpu_count

rule deepvariant__illumina_germline_variants:
    """
    Variant calling on Illumina reads using deep neural network.
    """
    input:
        bam   = f"{OUTDIR}/{{sample}}/minimap2/{{sample}}.illumina.sorted.bam",
        bai   = f"{OUTDIR}/{{sample}}/minimap2/{{sample}}.illumina.sorted.bam.bai",
        fasta = config["fasta"],
        fai   = config["fai"]
    output:
        vcf        = f"{OUTDIR}/{{sample}}/deepvariant/{{sample}}.illumina.deepvariant.vcf.gz",
        vcf_index  = f"{OUTDIR}/{{sample}}/deepvariant/{{sample}}.illumina.deepvariant.vcf.gz.tbi",
        gvcf       = f"{OUTDIR}/{{sample}}/deepvariant/{{sample}}.illumina.deepvariant.g.vcf.gz",
        gvcf_index = f"{OUTDIR}/{{sample}}/deepvariant/{{sample}}.illumina.deepvariant.g.vcf.gz.tbi",
        report     = f"{OUTDIR}/{{sample}}/deepvariant/{{sample}}.illumina.deepvariant.visual_report.html"
    params:
        model = "WGS"
    threads:
        lambda cores: cpu_count() - 2
    log:
        f"{OUTDIR}/{{sample}}/logs/deepvariant/{{sample}}.illumina.deepvariant.log",
    benchmark:
        f"{OUTDIR}/{{sample}}/logs/deepvariant/{{sample}}.illumina.deepvariant.benchmark"
    singularity:
        "docker://google/deepvariant:1.1.0"
    shell:
        """
        (/opt/deepvariant/bin/run_deepvariant \
            --model_type={params.model} \
            --ref={input.fasta} \
            --reads={input.bam} \
            --output_vcf={output.vcf} \
            --output_gvcf={output.gvcf} \
            --num_shards={threads}) > {log} 2>&1
        """


rule pepper__pacbio_germline_variants:
    """
    Variant calling on PacBio reads using deep neural network.
    """
    input:
        bam   = f"{OUTDIR}/{{sample}}/pbmm2/{{sample}}.pacbio.sorted.bam",
        bai   = f"{OUTDIR}/{{sample}}/pbmm2/{{sample}}.pacbio.sorted.bam.bai",
        fasta = config["fasta"],
        fai   = config["fai"]
    output:
        vcf           = f"{OUTDIR}/{{sample}}/pepper/{{sample}}.pacbio.pepper.vcf.gz",
        vcf_index     = f"{OUTDIR}/{{sample}}/pepper/{{sample}}.pacbio.pepper.vcf.gz.tbi",
        gvcf          = f"{OUTDIR}/{{sample}}/pepper/{{sample}}.pacbio.pepper.g.vcf.gz",
        gvcf_index    = f"{OUTDIR}/{{sample}}/pepper/{{sample}}.pacbio.pepper.g.vcf.gz.tbi",
        report        = f"{OUTDIR}/{{sample}}/deepvariant/{{sample}}.pacbio.pepper.visual_report.html"
    params:
        output_dir    = f"{OUTDIR}/{{sample}}/pepper",
        output_prefix = f"{{sample}}.pacbio.pepper"
    threads:
        lambda cores: cpu_count() - 2
    log:
        out = f"{OUTDIR}/{{sample}}/logs/pepper/{{sample}}.pacbio.pepper.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/logs/pepper/{{sample}}.pacbio.pepper.benchmark"
    singularity:
        "docker://kishwars/pepper_deepvariant:r0.4"
    shell:
        """
        (run_pepper_margin_deepvariant call_variant \
            -f {input.fasta} \
            -b {input.bam} \
            -o {params.output_dir} \
            -p {params.output_prefix} \
            -t {threads} \
            -s {wildcards.sample} \
            --gvcf \
            --ccs) > {log} 2>&1
        """


rule pepper__nanopore_germline_variants:
    """
    Variant calling using deep neural network.
    """
    input:
        bam   = f"{OUTDIR}/{{sample}}/minimap2/{{sample}}.nanopore.sorted.bam",
        bai   = f"{OUTDIR}/{{sample}}/minimap2/{{sample}}.nanopore.sorted.bam.bai",
        fasta = config["fasta"],
        fai   = config["fai"]
    output:
        vcf           = f"{OUTDIR}/{{sample}}/pepper/{{sample}}.nanopore.pepper.vcf.gz",
        vcf_index     = f"{OUTDIR}/{{sample}}/pepper/{{sample}}.nanopore.pepper.vcf.gz.tbi",
        gvcf          = f"{OUTDIR}/{{sample}}/pepper/{{sample}}.nanopore.pepper.g.vcf.gz",
        gvcf_index    = f"{OUTDIR}/{{sample}}/pepper/{{sample}}.nanopore.pepper.g.vcf.gz.tbi",
        report        = f"{OUTDIR}/{{sample}}/deepvariant/{{sample}}.nanopore.pepper.visual_report.html"
    params:
        output_dir    = f"{OUTDIR}/{{sample}}/pepper",
        output_prefix = f"{{sample}}.nanopore.pepper"
    threads:
        lambda cores: cpu_count() - 2
    log:
        out = f"{OUTDIR}/{{sample}}/logs/pepper/{{sample}}.nanopore.pepper.vcf.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/logs/pepper/{{sample}}.nanopore.pepper.vcf.benchmark"
    singularity:
        "docker://kishwars/pepper_deepvariant:r0.4"
    shell:
        """
        (run_pepper_margin_deepvariant call_variant \
            -f {input.fasta} \
            -b {input.bam} \
            -o {params.output_dir} \
            -p {params.output_prefix} \
            -t {threads} \
            -s {wildcards.sample} \
            --gvcf \
            --ont) > {log} 2>&1
        """
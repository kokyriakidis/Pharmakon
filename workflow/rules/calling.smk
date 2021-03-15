from multiprocessing import cpu_count

rule deepvariant__illumina_germline_variants:
    """
    Variant calling on Illumina reads using deep neural network.
    """
    input:
        bam   = f"{OUTDIR}/{{sample}}/minimap2/{{sample}}.{BUILD}.illumina.sorted.bam",
        bai   = f"{OUTDIR}/{{sample}}/minimap2/{{sample}}.{BUILD}.illumina.sorted.bam.bai",
        fasta = config["fasta"],
        fai   = config["fai"]
    output:
        vcf        = f"{OUTDIR}/{{sample}}/deepvariant/{{sample}}.{BUILD}.illumina.deepvariant.vcf.gz",
        vcf_index  = f"{OUTDIR}/{{sample}}/deepvariant/{{sample}}.{BUILD}.illumina.deepvariant.vcf.gz.tbi",
        gvcf       = f"{OUTDIR}/{{sample}}/deepvariant/{{sample}}.{BUILD}.illumina.deepvariant.g.vcf.gz",
        gvcf_index = f"{OUTDIR}/{{sample}}/deepvariant/{{sample}}.{BUILD}.illumina.deepvariant.g.vcf.gz.tbi",
        report     = f"{OUTDIR}/{{sample}}/deepvariant/{{sample}}.{BUILD}.illumina.deepvariant.visual_report.html"
    params:
        model = "WGS"
    threads:
        lambda cores: cpu_count() - 2
    log:
        f"{OUTDIR}/{{sample}}/deepvariant/logs/{{sample}}.{BUILD}.illumina.deepvariant.log",
    benchmark:
        f"{OUTDIR}/{{sample}}/deepvariant/logs/{{sample}}.{BUILD}.illumina.deepvariant.benchmark"
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


rule pepper__pacbio_haplotag_bam_step_1:
    """
    Variant calling on PacBio reads using deep neural network.
    Generation of haplotagged BAM file.
    """
    input:
        bam   = f"{OUTDIR}/{{sample}}/pbmm2/{{sample}}.{BUILD}.pacbio.sorted.bam",
        bai   = f"{OUTDIR}/{{sample}}/pbmm2/{{sample}}.{BUILD}.pacbio.sorted.bam.bai",
        fasta = config["fasta"],
        fai   = config["fai"]
    output:
        phased_bam_final       = f"{OUTDIR}/{{sample}}/pepper/{{sample}}.{BUILD}.pacbio.pepper.margin.phased.haplotagged.bam",
        phased_bam_index_final = f"{OUTDIR}/{{sample}}/pepper/{{sample}}.{BUILD}.pacbio.pepper.margin.phased.haplotagged.bam.bai"
    params:
        output_dir           = f"{OUTDIR}/{{sample}}/pepper",
        phased_bam           = f"{OUTDIR}/{{sample}}/pepper/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam",
        phased_bam_index     = f"{OUTDIR}/{{sample}}/pepper/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam.bai",
        intermediate_files   = f"{OUTDIR}/{{sample}}/pepper/intermediate_files",
        chuncks              = f"{OUTDIR}/{{sample}}/pepper/MARGIN_PHASED.PEPPER_SNP_MARGIN.chunks.csv",
        pepper_snp_vcf       = f"{OUTDIR}/{{sample}}/pepper/PEPPER_SNP_OUPUT.vcf.gz",
        pepper_snp_vcf_index = f"{OUTDIR}/{{sample}}/pepper/PEPPER_SNP_OUPUT.vcf.gz.tbi"
    threads:
        lambda cores: cpu_count() - 2
    benchmark:
        f"{OUTDIR}/{{sample}}/pepper/logs/{{sample}}.{BUILD}.pacbio.pepper.step1.benchmark"
    singularity:
        "docker://kishwars/pepper_deepvariant:r0.4"
    shell:
        """
        run_pepper_margin_deepvariant call_variant \
            -f {input.fasta} \
            -b {input.bam} \
            -o {params.output_dir} \
            -t {threads} \
            -s {wildcards.sample} \
            --ccs
        
        mv {params.phased_bam} {output.phased_bam_final}
        mv {params.phased_bam_index} {output.phased_bam_index_final}
        mv {params.chuncks} {params.intermediate_files}
        mv {params.pepper_snp_vcf} {params.intermediate_files}
        mv {params.pepper_snp_vcf_index} {params.intermediate_files}
        """

rule deepvariant__pacbio_germline_variants_step_2:
    """
    Variant calling on PacBio reads using deep neural network.
    """
    input:
        phased_bam   = f"{OUTDIR}/{{sample}}/pepper/{{sample}}.{BUILD}.pacbio.pepper.margin.phased.haplotagged.bam",
        phased_bai   = f"{OUTDIR}/{{sample}}/pepper/{{sample}}.{BUILD}.pacbio.pepper.margin.phased.haplotagged.bam.bai",
        fasta = config["fasta"],
        fai   = config["fai"]
    output:
        vcf           = f"{OUTDIR}/{{sample}}/pepper/{{sample}}.{BUILD}.pacbio.pepper.margin.deepvariant.vcf.gz",
        vcf_index     = f"{OUTDIR}/{{sample}}/pepper/{{sample}}.{BUILD}.pacbio.pepper.margin.deepvariant.vcf.gz.tbi",
        gvcf          = f"{OUTDIR}/{{sample}}/pepper/{{sample}}.{BUILD}.pacbio.pepper.margin.deepvariant.g.vcf.gz",
        gvcf_index    = f"{OUTDIR}/{{sample}}/pepper/{{sample}}.{BUILD}.pacbio.pepper.margin.deepvariant.g.vcf.gz.tbi",
        report        = f"{OUTDIR}/{{sample}}/pepper/{{sample}}.{BUILD}.pacbio.pepper.margin.deepvariant.visual_report.html"
    params:
        model = "PACBIO"
    threads:
        lambda cores: cpu_count() - 2
    log:
        f"{OUTDIR}/{{sample}}/pepper/logs/3_deepvariant.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/pepper/logs/{{sample}}.{BUILD}.pacbio.pepper.step2.benchmark"
    singularity:
        "docker://google/deepvariant:1.1.0"
    shell:
        """
        (/opt/deepvariant/bin/run_deepvariant \
            --model_type={params.model} \
            --ref={input.fasta} \
            --reads={input.phased_bam} \
            --output_vcf={output.vcf} \
            --output_gvcf={output.gvcf} \
            --num_shards={threads} \
            --use_hp_information) > {log} 2>&1
        """

rule margin__pacbio_phase_vcf_step_3:
    """
    VCF phasing using Margin.
    """
    input:
        phased_bam             = f"{OUTDIR}/{{sample}}/pepper/{{sample}}.{BUILD}.pacbio.pepper.margin.phased.haplotagged.bam",
        phased_bam_index       = f"{OUTDIR}/{{sample}}/pepper/{{sample}}.{BUILD}.pacbio.pepper.margin.phased.haplotagged.bam.bai",
        fasta                  = config["fasta"],
        fai                    = config["fai"],
        vcf                    = f"{OUTDIR}/{{sample}}/pepper/{{sample}}.{BUILD}.pacbio.pepper.margin.deepvariant.vcf.gz",
        vcf_index              = f"{OUTDIR}/{{sample}}/pepper/{{sample}}.{BUILD}.pacbio.pepper.margin.deepvariant.vcf.gz.tbi",
        gvcf                   = f"{OUTDIR}/{{sample}}/pepper/{{sample}}.{BUILD}.pacbio.pepper.margin.deepvariant.g.vcf.gz",
        gvcf_index             = f"{OUTDIR}/{{sample}}/pepper/{{sample}}.{BUILD}.pacbio.pepper.margin.deepvariant.g.vcf.gz.tbi",
    output:
        phased_vcf       = f"{OUTDIR}/{{sample}}/pepper/{{sample}}.{BUILD}.pacbio.pepper.margin.deepvariant.phased.vcf.gz",
        phased_vcf_index = f"{OUTDIR}/{{sample}}/pepper/{{sample}}.{BUILD}.pacbio.pepper.margin.deepvariant.phased.vcf.gz.tbi",
        phaseset_bed     = f"{OUTDIR}/{{sample}}/pepper/{{sample}}.{BUILD}.pacbio.pepper.margin.deepvariant.phaseset.bed"
    params:
        output_prefix      = f"{OUTDIR}/{{sample}}/pepper/{{sample}}.{BUILD}.pacbio.pepper.margin.deepvariant",
        intermediate_files = f"{OUTDIR}/{{sample}}/pepper/intermediate_files"
    threads:
        lambda cores: cpu_count() - 2
    benchmark:
        f"{OUTDIR}/{{sample}}/pepper/logs/{{sample}}.{BUILD}.pacbio.pepper.step3.benchmark"
    singularity:
        "docker://kishwars/pepper_deepvariant:r0.4"
    shell:
        """
        margin phase \
            {input.phased_bam} \
            {input.fasta} \
            {input.vcf} \
            /opt/margin_dir/params/misc/allParams.phase_vcf.json \
            -t {threads} \
            -M \
            -o {params.output_prefix}

        rm  {params.output_prefix}.chunks.csv
        bgzip -c {params.output_prefix}.phased.vcf > {output.phased_vcf}
        rm {params.output_prefix}.phased.vcf
        tabix -p vcf {output.phased_vcf}
        rm -rf {params.intermediate_files}
        """


rule pepper__nanopore_germline_variants:
    """
    Variant calling on Nanopore reads using deep neural network.
    """
    input:
        bam   = f"{OUTDIR}/{{sample}}/minimap2/{{sample}}.{BUILD}.nanopore.sorted.bam",
        bai   = f"{OUTDIR}/{{sample}}/minimap2/{{sample}}.{BUILD}.nanopore.sorted.bam.bai",
        fasta = config["fasta"],
        fai   = config["fai"]
    output:
        vcf                = f"{OUTDIR}/{{sample}}/pepper/{{sample}}.{BUILD}.nanopore.pepper.margin.deepvariant.vcf.gz",
        vcf_index          = f"{OUTDIR}/{{sample}}/pepper/{{sample}}.{BUILD}.nanopore.pepper.margin.deepvariant.vcf.gz.tbi",
        gvcf               = f"{OUTDIR}/{{sample}}/pepper/{{sample}}.{BUILD}.nanopore.pepper.margin.deepvariant.g.vcf.gz",
        gvcf_index         = f"{OUTDIR}/{{sample}}/pepper/{{sample}}.{BUILD}.nanopore.pepper.margin.deepvariant.g.vcf.gz.tbi",
        report             = f"{OUTDIR}/{{sample}}/pepper/{{sample}}.{BUILD}.nanopore.pepper.margin.deepvariant.visual_report.html",
        phased_vcf         = f"{OUTDIR}/{{sample}}/pepper/{{sample}}.{BUILD}.nanopore.pepper.margin.deepvariant.phased.vcf.gz",
        phased_vcf_index   = f"{OUTDIR}/{{sample}}/pepper/{{sample}}.{BUILD}.nanopore.pepper.margin.deepvariant.phased.vcf.gz.tbi",
        phased_bam         = f"{OUTDIR}/{{sample}}/pepper/{{sample}}.{BUILD}.nanopore.pepper.margin.phased.haplotagged.bam",
        phased_bam_index   = f"{OUTDIR}/{{sample}}/pepper/{{sample}}.{BUILD}.nanopore.pepper.margin.phased.haplotagged.bam.bai",
        phaseset_bed       = f"{OUTDIR}/{{sample}}/pepper/{{sample}}.{BUILD}.nanopore.pepper.margin.deepvariant.phaseset.bed"
    params:
        output_dir    = f"{OUTDIR}/{{sample}}/pepper",
        output_prefix = f"{{sample}}.{BUILD}.nanopore.pepper.margin.deepvariant",
        intermediate_files = f"{OUTDIR}/{{sample}}/pepper/intermediate_files"
    threads:
        lambda cores: cpu_count() - 2
    benchmark:
        f"{OUTDIR}/{{sample}}/pepper/logs/{{sample}}.{BUILD}.nanopore.pepper.margin.deepvariant.benchmark"
    singularity:
        "docker://kishwars/pepper_deepvariant:r0.4"
    shell:
        """
        run_pepper_margin_deepvariant call_variant \
            -b {input.bam} \
            -f {input.fasta} \
            -o {params.output_dir} \
            -p {params.output_prefix} \
            -t {threads} \
            -s {wildcards.sample} \
            --gvcf \
            --phased_output \
            --ont

        mv {params.intermediate_files}/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam {output.phased_bam}
        mv {params.intermediate_files}/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam.bai {output.phased_bam_index}
        rm -rf {params.intermediate_files}
        """


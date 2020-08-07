from sys import platform as _platform
genome_build = config["ref"]["build"]

rule get_genome:
    output:
        "reference/{genome_build}/{genome_build}.fasta".format(genome_build=genome_build)
    log:
        "logs/reference/get-{genome_build}.log".format(genome_build=genome_build)
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"]
    cache: True
    wrapper:
        "0.64.0/bio/reference/ensembl-sequence"


checkpoint genome_faidx:
    input:
        "reference/{genome_build}/{genome_build}.fasta".format(genome_build=genome_build)
    output:
        "reference/{genome_build}/{genome_build}.fasta.fai".format(genome_build=genome_build)
    log:
        "logs/reference/{genome_build}-faidx.log".format(genome_build=genome_build)
    cache: True
    wrapper:
        "0.64.0/bio/samtools/faidx"


rule genome_dict:
    input:
        "reference/{genome_build}/{genome_build}.fasta".format(genome_build=genome_build)
    output:
        "reference/{genome_build}/{genome_build}.dict".format(genome_build=genome_build)
    log:
        "logs/reference/create_dict.log"
    conda:
        "../envs/samtools.yaml"
    cache: True
    shell:
        "samtools dict {input} > {output} 2> {log} "


if _platform == "darwin":

    rule get_known_variation:
        input:
            # use fai to annotate contig lengths for GATK BQSR
            fai=rules.genome_faidx.output[0]
        output:
            vcf="reference/{genome_build}/{genome_build}_variation.vcf.gz".format(genome_build=genome_build)
        log:
            "logs/reference/get-known-variants.log"
        params:
            species=config["ref"]["species"],
            build=config["ref"]["build"],
            release=config["ref"]["release"],
            type="all"
        cache: True
        wrapper:
            "0.64.0/bio/reference/ensembl-variation"


    rule remove_iupac_codes:
        input:
            rules.get_known_variation.output["vcf"]
        output:
            "reference/{genome_build}/{genome_build}_variation_noiupac.vcf.gz".format(genome_build=genome_build)
        log:
            "logs/reference/fix-iupac-alleles.log"
        conda:
            "../envs/rbt.yaml"
        cache: True
        shell:
            "rbt vcf-fix-iupac-alleles < {input} | bcftools view -Oz > {output}"


    rule tabix_known_variants:
        input:
            rules.remove_iupac_codes.output[0]
        output:
            "reference/{genome_build}/{genome_build}_variation_noiupac.vcf.gz.tbi".format(genome_build=genome_build)
        log:
            "logs/reference/tabix/variation.log"
        params:
            "-p vcf"
        cache: True
        wrapper:
            "0.64.0/bio/tabix"

    rule bwa_index:
        input:
            "reference/{genome_build}/{genome_build}.fasta".format(genome_build=genome_build)
        output:
            "reference/{genome_build}/{genome_build}.amb".format(genome_build=genome_build),
            "reference/{genome_build}/{genome_build}.ann".format(genome_build=genome_build),
            "reference/{genome_build}/{genome_build}.bwt".format(genome_build=genome_build),
            "reference/{genome_build}/{genome_build}.pac".format(genome_build=genome_build),
            "reference/{genome_build}/{genome_build}.sa".format(genome_build=genome_build)
        log:
            "logs/reference/bwa_index.log"
        resources:
            mem_mb=369000
        cache: True
        wrapper:
            "0.64.0/bio/bwa/index"

elif _platform == "linux" or _platform == "linux2":

    rule bwa_mem2_index:
        input:
            "reference/{genome_build}/{genome_build}.fasta".format(genome_build=genome_build)
        output:
            "reference/{genome_build}/{genome_build}.0123".format(genome_build=genome_build),
            "reference/{genome_build}/{genome_build}.amb".format(genome_build=genome_build),
            "reference/{genome_build}/{genome_build}.ann".format(genome_build=genome_build),
            "reference/{genome_build}/{genome_build}.bwt.2bit.64".format(genome_build=genome_build),
            "reference/{genome_build}/{genome_build}.bwt.8bit.32".format(genome_build=genome_build),
            "reference/{genome_build}/{genome_build}.pac".format(genome_build=genome_build)
        log:
            "logs/reference/bwa_mem2_index.log"
        resources:
            mem_mb=369000
        cache: True
        wrapper:
            "0.64.0/bio/bwa-mem2/index"


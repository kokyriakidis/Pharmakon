from sys import platform as _platform
genome_build = config["ref"]["build"]


rule get_genome:
    output:
        "resources/{genome_build}/{genome_build}.fasta".format(genome_build=genome_build)
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
        "resources/{genome_build}/{genome_build}.fasta".format(genome_build=genome_build)
    output:
        "resources/{genome_build}/{genome_build}.fasta.fai".format(genome_build=genome_build)
    log:
        "logs/reference/{genome_build}-faidx.log".format(genome_build=genome_build)
    cache: True
    wrapper:
        "0.64.0/bio/samtools/faidx"


rule genome_dict:
    input:
        "resources/{genome_build}/{genome_build}.fasta".format(genome_build=genome_build)
    output:
        "resources/{genome_build}/{genome_build}.dict".format(genome_build=genome_build)
    log:
        "logs/reference/create_dict.log"
    conda:
        "../envs/samtools.yaml"
    cache: True
    shell:
        "samtools dict {input} > {output} 2> {log} "


if _platform == "darwin":
    rule bwa_index:
        input:
            "resources/{genome_build}/{genome_build}.fasta".format(genome_build=genome_build)
        output:
            "resources/{genome_build}/{genome_build}.amb".format(genome_build=genome_build),
            "resources/{genome_build}/{genome_build}.ann".format(genome_build=genome_build),
            "resources/{genome_build}/{genome_build}.bwt".format(genome_build=genome_build),
            "resources/{genome_build}/{genome_build}.pac".format(genome_build=genome_build),
            "resources/{genome_build}/{genome_build}.sa".format(genome_build=genome_build)
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
            "resources/{genome_build}/{genome_build}.fasta".format(genome_build=genome_build)
        output:
            "resources/{genome_build}/{genome_build}.0123".format(genome_build=genome_build),
            "resources/{genome_build}/{genome_build}.amb".format(genome_build=genome_build),
            "resources/{genome_build}/{genome_build}.ann".format(genome_build=genome_build),
            "resources/{genome_build}/{genome_build}.bwt.2bit.64".format(genome_build=genome_build),
            "resources/{genome_build}/{genome_build}.bwt.8bit.32".format(genome_build=genome_build),
            "resources/{genome_build}/{genome_build}.pac".format(genome_build=genome_build)
        log:
            "logs/reference/bwa_mem2_index.log"
        resources:
            mem_mb=369000
        cache: True
        wrapper:
            "0.64.0/bio/bwa-mem2/index"


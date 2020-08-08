genome_build = config["ref"]["build"]

rule get_genome2:
    output:
        "reference/{genome_build}/{genome_build}.fa".format(genome_build=genome_build),
        #"reference/{genome_build}/{genome_build}.fa.fai".format(genome_build=genome_build),
        #"reference/{genome_build}/{genome_build}.fa.gz".format(genome_build=genome_build),
        #"reference/{genome_build}/{genome_build}.fa.gz.fai".format(genome_build=genome_build),
        #"reference/{genome_build}/{genome_build}.fa.gz.gzi".format(genome_build=genome_build),
        #"reference/{genome_build}/{genome_build}.dict".format(genome_build=genome_build)
    log:
        "logs/reference/get-{genome_build}.log".format(genome_build=genome_build)
    params:
        genome_build = config["ref"]["build"]
    shell:
        "mkdir reference/{params.genome_build}/"



        #wget --no-check-certificate -c https://s3.amazonaws.com/biodata/genomes/hg19-seq.tar.gz
        #tar -xzvpf hg19-seq.tar.gz
        #gunzip -c seq/hg19.fa.gz > seq/hg19.fa
        #touch seq/hg19.fa.fai
        #touch seq/hg19.dict
        #""""

        #"(bwa mem -R '{params.rg}' -t {threads} {input} | "
        #"samtools view -Sb - > {output}) 2> {log}"
def get_regions_param(regions=config["processing"].get("restrict-regions"), default=""):
    if regions:
        params = "--regions '{}' ".format(regions)
        return params
    return default

rule deepvariant:
    input:
        bam="dedup/{sample}/{sample}.bam",
        ref=rules.get_genome.output[0]
    output:
        vcf="calls/{sample}/{sample}.vcf.gz"
    params:
        model=config["params"]["deepvariant"]["model"],
        extra=get_regions_param
    threads: 2
    log:
        "logs/deepvariant/{sample}/{sample}/stdout.log"
    wrapper:
        "0.64.0/bio/deepvariant"



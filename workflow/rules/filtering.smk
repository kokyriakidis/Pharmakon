if _platform == "darwin":

    def get_vartype_arg(wildcards):
        return "--select-type-to-include {}".format(
            "SNP" if wildcards.vartype == "snvs" else "INDEL")


    rule select_calls:
        input:
            ref=rules.get_genome.output[0],
            vcf=rules.merge_variants.output[0]
        output:
            vcf=temp("{project_dir}/{sample}/filtered/all.{vartype}.vcf.gz")
        params:
            extra=get_vartype_arg
        log:
            "{project_dir}/{sample}/logs/gatk/selectvariants/{vartype}.log"
        wrapper:
            "0.59.0/bio/gatk/selectvariants"


    def get_filter(wildcards):
        return {
            "snv-hard-filter":
            config["filtering"]["hard"][wildcards.vartype]}


    rule hard_filter_calls:
        input:
            ref=rules.get_genome.output[0],
            vcf=rules.select_calls.output[0]
        output:
            vcf=temp("{project_dir}/{sample}/filtered/all.{vartype}.hardfiltered.vcf.gz")
        params:
            filters=get_filter
        log:
            "{project_dir}/{sample}/logs/gatk/variantfiltration/{vartype}.log"
        wrapper:
            "0.59.2/bio/gatk/variantfiltration"


    rule recalibrate_calls:
        input:
            vcf=rules.select_calls.output[0]
        output:
            vcf=temp("{project_dir}/{sample}/filtered/all.{vartype}.recalibrated.vcf.gz")
        params:
            extra=config["params"]["gatk"]["VariantRecalibrator"]
        log:
            "{project_dir}/{sample}/logs/gatk/variantrecalibrator/{vartype}.log"
        wrapper:
            "0.59.2/bio/gatk/variantrecalibrator"


    rule merge_calls:
        input:
            vcfs=expand("{project_dir}/{sample}/filtered/all.{vartype}.{filtertype}.vcf.gz",
                    vartype=["snvs", "indels"],
                    filtertype="recalibrated"
                                if config["filtering"]["vqsr"]
                                else "hardfiltered")
        output:
            vcf="{project_dir}/{sample}/final_vcf/{sample}.vcf.gz"
        log:
            "{project_dir}/{sample}/logs/picard/merge-filtered.log"
        wrapper:
            "0.59.2/bio/picard/mergevcfs"
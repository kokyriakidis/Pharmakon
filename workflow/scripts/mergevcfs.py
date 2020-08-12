##############################
from typing import Dict, List, Optional
import pysam
import pandas as pd
import numpy as np
from snakemake.shell import shell


def read_gene_table(fn: str) -> Dict[str, Dict[str, str]]:
    """
    Read gene table file.

    Returns:
        dict[str, dict[str, str]]: Gene table object.

    Args:
        fn (str): Gene table file.

    Examples:

        >>> gene_table = read_gene_table("gene_table.txt")
        >>> print(gene_table["cyp2d6"]["hg19_start"])
        42522500
    """

    result = {}

    with open(fn) as f:
        header = next(f).strip().split("\t")

        for line in f:
            fields = line.strip().split("\t")
            name = fields[1]
            result[name] = dict(zip(header, fields))

    return result

def get_gene_table() -> Dict[str, Dict[str, str]]:
    """
    Get gene table object.

    Returns:
        dict[str, dict[str, str]]: Gene table object.
    """

    return read_gene_table("workflow/resources/stargazer/gene_table.txt")

def get_target_genes() -> List[str]:
    """Get the list of target gene names.

    Returns:
        list[str]: A list of gene names.
    """
    gene_table = get_gene_table()
    return [k for k, v in gene_table.items() if v["type"] == "target"]

def get_target_region(tg: str, gb: str) -> str:
    """Get the genomic region for the target gene.

    Returns:
        str: Genomic region.

    Args:
        tg (str): Target gene.
        gb (str): Genome build (hg19, hg38).
    """
    gene_table = get_gene_table()
    target_genes = [k for k, v in gene_table.items() if v["type"] == "target"]

    if tg not in target_genes:
        raise ValueError(f"'{tg}' is not among target genes: {target_genes}")

    return gene_table[tg][f"{gb}_region"]


def is_chr(bam: str) -> bool:
    """
    Check whether SN tags in BAM file contain "chr" string.

    Returns:
        bool: True if found.

    Args:
        bam (str): BAM file.
    """

    header = pysam.view("-H", bam).strip().split("\n")

    l = []

    for line in header:
        fields = line.split("\t")
        if "@SQ" == fields[0]:
            for field in fields:
                if "SN:" in field:
                    l.append(field.replace("SN:", ""))

    return any(["chr" in x for x in l])

##############################


##### GET SELECTED GENES #####
stargazer_target_genes = get_target_genes()
    
if snakemake.config["params"]["stargazer"]["target_genes"] == "ALL":
    selected_genes = stargazer_target_genes
else:
    selected_genes = []
    for gene in snakemake.config["params"]["stargazer"]["target_genes"].split(","):
        selected_genes.append(gene.strip().lower())
    for gene in selected_genes:
        if gene not in stargazer_target_genes:
            raise ValueError(f"Unrecognized target gene found: {gene}")


##### GET PROJECT DIR #####
project_dir=snakemake.config["project_dir"]

##### GET SAMPLES #####
samples = pd.read_table(snakemake.config["samples"]).set_index("sample_name", drop=False)

##### GET FILTERTYPE #####
filtertype="recalibrated" if snakemake.config["filtering"]["vqsr"] else "hardfiltered"


for sample in samples["sample_name"]:
    for gene in selected_genes:
        INPUT1 = project_dir + "/" + sample + "/genes/" + gene + "/filtered_vcf/" + gene + ".snvs." + filtertype + ".vcf.gz"
        INPUT2 = project_dir + "/" + sample + "/genes/" + gene + "/filtered_vcf/" + gene + ".indels." + filtertype + ".vcf.gz"
        output = project_dir + "/" + sample + "/genes/" + gene + "/final_vcf/" + gene + ".vcf.gz"
        log = snakemake.log_fmt_shell(stdout=False, stderr=True)

        shell("picard MergeVcfs I={INPUT1} I={INPUT2} OUTPUT={output} {log}")
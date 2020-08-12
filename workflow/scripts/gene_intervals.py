##############################
from sys import platform as _platform
from typing import Dict, List, Optional
import pysam
import pandas as pd
import numpy as np

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

if snakemake.config["params"]["stargazer"]["target_genes"] == "ALL":
    selected_genes = stargazer_target_genes
else:
    selected_genes = []
    for gene in config["params"]["stargazer"]["target_genes"].split(","):
        selected_genes.append(gene.strip().lower())
    for gene in selected_genes:
        if gene not in stargazer_target_genes:
            raise ValueError(f"Unrecognized target gene found: {gene}")

project_dir=snakemake.config["project_dir"]
samples = pd.read_table(snakemake.config["samples"]).set_index("sample_name", drop=False)
gene_table = get_gene_table()
if snakemake.config["ref"]["build"] == "GRCh38":
    genome_build = "hg38"
elif snakemake.config["ref"]["build"] == "GRCh37":
    genome_build = "hg19"
for sample in samples["sample_name"]:
    if _platform == "darwin" or config["variant_tool"] == "gatk":
        # case 1: no duplicate removal
        bam = f"{project_dir}/{sample}/applybqsr/{sample}.bam"
    else:
        # case 1: no duplicate removal
        bam = f"{project_dir}/{sample}/mapped/{sample}.sorted.bam"
        if snakemake.config["processing"]["remove-duplicates"]:
            # case 2: remove duplicates
            bam = f"{project_dir}/{sample}/dedup/{sample}.bam"
            
    _ = is_chr(bam)

    if _:
        chr_str = "chr"
    elif not _:
        chr_str = ""
    else:
        raise ValueError("Mixed types of SN tags found.") 

    target_region = ""
    for gene in selected_genes:
        target_region += chr_str + "{chromosome}\t{start}\t{end}\n".format(chromosome=gene_table[gene]["chr"].replace("chr", ""), start=gene_table[gene]["{genome_build}_start".format(genome_build=genome_build)], end=gene_table[gene]["{genome_build}_end".format(genome_build=genome_build)])
        with open(
            f"{project_dir}/{sample}/intervals/{gene}.bed", "w"
        ) as f:
            f.write(target_region)
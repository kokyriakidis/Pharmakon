import os
import pysam
import statistics
from typing import Optional, List, TextIO

def get_gene_table() -> Dict[str, Dict[str, str]]:
    """
    Get gene table object.

    Returns:
        dict[str, dict[str, str]]: Gene table object.
    """

    return read_gene_table("workflow/resources/stargazer/gene_table.txt")


def read_gene_table(
        fn: str
    ) -> Dict[str, Dict[str, str]]:
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


def sm_tag(bam: str) -> str:
    """
    Extract SM tag from BAM file.

    Returns:
        str: SM tag.

    Args:
        bam (str): BAM file.
    """

    header = pysam.view("-H", bam).strip().split("\n")

    l = []

    for line in header:
        fields = line.split("\t")
        if "@RG" == fields[0]:
            for field in fields:
                if "SM:" in field:
                    l.append(field.replace("SM:", ""))

    l = list(set(l))

    if not l:
        raise ValueError(f"SM tag not found: {bam}")

    if len(l) > 1:
        result = l[0]
    else:
        result = l[0]

    return result


def parse_region(
        region: str,
        omit: bool = False
    ) -> List[str]:
    """
    Parse region.

    Returns:
        list[str]: Parsed region [chr, start, end].

    Args:
        region (str): Region to be parsed.
        omit (bool): Remove the 'chr' string.
    """
    if omit:
        chr = region.split(":")[0].replace("chr", "")
    else:
        chr = region.split(":")[0]

    return (
        chr,
        int(region.split(":")[1].split("-")[0]),
        int(region.split(":")[1].split("-")[1]),
    )


def sort_regions(regions: List[str]) -> List[str]:
    """
    Sort regions.

    Returns:
        list[str]: Sorted regions.

    Args:
        regions (list[str]): Regions.
    """

    def f(x):
        r = parse_region(x)
        if "X" in r[0]:
            chr = 23
        elif "Y" in r[0]:
            chr = 24
        else:
            chr = int(r[0].replace("chr", ""))
        return (chr, r[1], r[2])
    return sorted(regions, key = f)


def bam2sdf(
        gb:str,
        tg: str,
        cg: str,
        bam: List[str],
        **kwargs
    ) -> str:
    """
    Create SDF file from BAM file(s).

    Returns:
        str: SDF file.

    Args:
        gb (str): Genome build (hg19, hg38).
        tg (str): Target gene.
        cg (str): Control gene or region.
        bam (list[str]): BAM file(s).
    """

    gene_table = get_gene_table()

    targets = [k for k, v in gene_table.items() if v["type"] == "target"]

    if tg not in targets:
        raise ValueError(f"'{tg}' is not among target genes: {targets}")

    tr = gene_table[tg][f"{gb}_region"].replace("chr", "")

    if "chr" in cg or ":" in cg:
        cr = cg.replace("chr", "")

    else:
        controls = [k for k, v in gene_table.items() if v["control"] == "yes"]

        if cg not in controls:
            raise ValueError(f"'{cg}' is not among control genes: {controls}")

        cr = gene_table[cg][f"{gb}_region"].replace("chr", "")

    regions = sort_regions([tr, cr])

    # Get sample and sequence names from BAM headers.
    sm = []
    sn = []
    for x in bam:
        sm.append(sm_tag(x))

        result = pysam.view("-H", x).strip().split("\n")
        for line in result:
            fields = line.split("\t")
            if "@SQ" == fields[0]:
                for field in fields:
                    if "SN:" in field:
                        y = field.replace("SN:", "")
                        if y not in sn:
                            sn.append(y)

    # Determine whether the "chr" string should be used.
    if any(["chr" in x for x in sn]):
        chr_str = "chr"
    else:
        chr_str = ""

    result = ""

    for region in regions:
        temp = pysam.depth("-a", "-Q", "1", "-r", f"{chr_str}{region}", *bam)
        result += temp

    return result


def sdf2gdf(
        fn: str,
        id: List[str],
        f: Optional[TextIO] = None, 
        **kwargs
    ) -> str:
    """
    Create GDF file from SDF file.

    Returns:
        str: GDF file.

    Args:
        fn (str): SDF file.
        id (list[str]): Sample ID(s).
        f (TextIO, optional): SDF file.
    """

    if fn:
        f = open(fn)

    # Get the header.
    result = "Locus\tTotal_Depth\tAverage_Depth_sample"
    for x in id:
        result += f"\tDepth_for_{x}"
    result += "\n"

    # Check the sample count with the first line.
    fields1 = next(f).strip().split("\t")
    if len(fields1) - 2 != len(id):
        raise ValueError("incorrect sample count")
    locus = f"{fields1[0]}:{fields1[1]}"
    depth = [int(x) for x in fields1[2:]]
    avg = round(statistics.mean(depth), 2)
    fields2 = [locus, sum(depth), avg] + depth
    result += "\t".join([str(x) for x in fields2]) + "\n"

    # Read rest of the lines.
    for line in f:
        fields1 = line.strip().split("\t")
        locus = f"{fields1[0]}:{fields1[1]}"
        depth = [int(x) for x in fields1[2:]]
        avg = round(statistics.mean(depth), 2)
        fields2 = [locus, sum(depth), avg] + depth
        result += "\t".join([str(x) for x in fields2]) + "\n"

    if fn:
        f.close()

    return result



def bam2gdf(
        genome_build: str,
        target_gene : str,
        control_gene: str,
        output_file: str,
        bam_file: List[str],
        bam_dir: Optional[str] = None,
        bam_list: Optional[str] = None,
        **kwargs
    ) -> None:
    """Convert BAM files to a GDF file.

    This command calculates read depth from BAM files and then outputs a
    GDF (GATK-DepthOfCoverage Format) file, which is one of the input 
    files for the Stargazer program. Even though ``gatk DepthOfCoverage`` 
    could still be used to make GDF files, we recommend that you use this 
    command because the former is too heavy (i.e. requires too much memory) 
    for such a simple task (i.e. counting reads). The latter uses 
    ``samtools depth`` under the hood, which is way faster and requires 
    way less memory. Another nice about using ``bam2gdf`` instead of 
    ``samtools depth`` is that everything is already parametrized for 
    compatibility with Stargazer. 

    .. note::
        You do NOT need to install ``samtools`` to run this command.

    Args:
        genome_build (str):
            Genome build ('hg19' or 'hg38').
        target_gene (str):
            Name of target gene (e.g. 'cyp2d6').
        control_gene (str):
            Name or region of control gene (e.g. ‘vdr’, 
            ‘chr12:48232319-48301814’)
        output_file (str):
            Write output to this file.
        bam_file (list[str]):
            Input BAM files.
        bam_dir (str, optional):
            Use all BAM files in this directory as input.
        bam_list (str, optional):
            List of input BAM files, one file per line.
    """
    # Parse keyward arguments from the decorator.
    input_files = kwargs["input_files"]

    sdf = bam2sdf(genome_build, target_gene, control_gene, input_files)
    sm = [sm_tag(x) for x in input_files]
    result = sdf2gdf(sdf, sm)
    with open(output_file, "w") as f:
        f.write(result)


bam2gdf(
        snakemake.input["genome_build"],
        snakemake.input["target_gene"],
        snakemake.input["control_gene"],
        snakemake.output["output_file"],
        bam_file = snakemake.input["bam_file"]
    )
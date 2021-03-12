import os
import csv
import sys
import datetime
import pkgutil
from typing import TextIO, Optional

from .common import get_target_genes
from .gt2pt import phenotyper

def _add_overview_section(genotype_table, pair_table, target_genes):
    table = (
        "<table>\n"
        "<tr>\n"
        "  <th style='width: 10%;'>No.</th>\n"
        "  <th style='width: 10%;'>Gene</th>\n"
        "  <th style='width: 20%;'>Genotype</th>\n"
        "  <th style='width: 10%;'>Total AS</th>\n"
        "  <th style='width: 30%;'>Predicted Phenotype</th>\n"
        "  <th style='width: 10%;'>Drugs</th>\n"
        "  <th style='width: 10%;'>Guidelines</th>\n"
        "</tr>\n"
    )
    words = ["normal", "unknown", "-"]
    for i, gene in enumerate(target_genes, 1):
        if any([x in genotype_table[gene]["phenotype"] for x in words]):
            color = "black"
        else:
            color = "red"
        genotype = genotype_table[gene]["hap1_main"] + "/" + genotype_table[gene]["hap2_main"]
        score = genotype_table[gene]["dip_score"]
        phenotype = genotype_table[gene]["phenotype"]
        drugs = len([x["Drug"] for x in pair_table if x["Gene"] == gene.upper()])
        guidelines = len([x for x in pair_table if x["Gene"] == gene.upper() and x["Guideline"] != "-"])
        table += (
            f"<tr style='color: {color};'>\n"
            f"  <td>{i}</td>\n"
            f"  <td>{gene.upper()}</td>\n"
            f"  <td>{genotype}</td>\n"
            f"  <td>{score}</td>\n"
            f"  <td>{phenotype}</td>\n"
            f"  <td>{drugs}</td>\n"
            f"  <td>{guidelines}</td>\n"
            "</tr>\n"
        )
    return table + "</table>"

def _add_genotypes_section(genotype_table, target_genes):
    table = (
        "<table>\n"
        "<tr>\n"
        "  <th style='width: 5%;'>No.</th>\n"
        "  <th style='width: 5%;'>Gene</th>\n"
        "  <th style='width: 20%;'>Star Allele</th>\n"
        "  <th style='width: 10%;'>AS</th>\n"
        "  <th style='width: 50%;'>SNVs/Indels</th>\n"
        "  <th style='width: 10%;'>SVs</th>\n"
        "</tr>\n"
    )
    for i, gene in enumerate(target_genes, 1):
        hmc1 = "<br />".join(genotype_table[gene]["hap1_main_core"].split(","))
        hmc2 = "<br />".join(genotype_table[gene]["hap2_main_core"].split(","))
        hap1_main = genotype_table[gene]["hap1_main"]
        hap1_score = genotype_table[gene]["hap1_score"]
        hap1_sv = genotype_table[gene]["hap1_sv"]
        hap2_main = genotype_table[gene]["hap2_main"]
        hap2_score = genotype_table[gene]["hap2_score"]
        hap2_sv = genotype_table[gene]["hap2_sv"]
        table += (
            "<tr>\n"
            f"  <td rowspan='2'>{i}</td>\n"
            f"  <td rowspan='2'>{gene.upper()}</td>\n"
            f"  <td>{hap1_main}</td>\n"
            f"  <td>{hap1_score}</td>\n"
            f"  <td>{hmc1}</td>\n"
            f"  <td>{hap1_sv}</td>\n"
            "</tr>\n"
            "<tr>\n"
            f"  <td>{hap2_main}</td>\n"
            f"  <td>{hap2_score}</td>\n"
            f"  <td>{hmc2}</td>\n"
            f"  <td>{hap2_sv}</td>\n"
            "</tr>\n"
        )
    return table + "</table>"

def _add_drugs_section(genotype_table, pair_table):
    table = (
        "<table>\n"
        "<tr>\n"
        "  <th style='width: 5%;'>No.</th>\n"
        "  <th style='width: 15%;'>Drug</th>\n"
        "  <th style='width: 5%;'>Gene</th>\n"
        "  <th style='width: 5%;'>Level</th>\n"
        "  <th style='width: 10%;'>FDA</th>\n"
        "  <th style='width: 60%;'>Guideline</th>\n"
        "</tr>\n"
    )
    for i, pair in enumerate(sorted(
        pair_table, key = lambda x: (x["Drug"].lower(), x["Gene"])), 1):

        gn = pair["Gene"].lower() # gene name

        # Gene-drug pairs are bolded if a genotype is available for the gene.
        if gn in genotype_table and genotype_table[gn]["hap1_main"] != "-":
            bold = "bold"
        else:
            bold = "normal"

        # Gene-drug pairs are in red if altered phenotype is predicted.
        if bold == "normal":
            color = "black"
        elif any([x in genotype_table[gn]["phenotype"]
            for x in ["normal", "unknown", "-"]]):
            color = "black"
        else:
            color = "red"

        table += (
        f"<tr style='font-weight: {bold}; color: {color};'>\n"
        f"  <td>{i}</td>\n"
        f"  <td>{pair['Drug']}</td>\n"
        f"  <td>{pair['Gene']}</td>\n"
        f"  <td>{pair['CPIC Level']}</td>\n"
        f"  <td>{pair['PGx on FDA Label']}</td>\n"
        f"  <td>{pair['Guideline']}</td>\n"
        "</tr>\n"
        )
    return table + "</table>"

def _add_recommendations_section(action_table):
    string = ""
    for chemical in action_table:
        for gene in action_table[chemical]:
            description = (
                f"{action_table[chemical][gene]['summary']} "
                f"[PharmGKB Link: {action_table[chemical][gene]['url']}]"
            )
            table = (
                f"<h3>{chemical}-{gene}</h3>\n"
                f"<p>{description}</p>\n"
                "<table>\n"
                "<tr>\n"
                "  <th style='width: 20%;'>Phenotype</th>\n"
                "  <th style='width: 80%;'>Recommendation</th>\n"
                "</tr>\n"
            )
            for phenotype in action_table[chemical][gene]["pt"]:
                table += (
                    "<tr>\n"
                    f"  <td>{phenotype}</td>\n"
                    f"  <td>{action_table[chemical][gene]['pt'][phenotype]}</td>\n"
                    "</tr>\n"
                )
            table += "</table>\n"
            string += table
    return string

def _read_action_table():
    result = {}
    p = os.path.dirname(__file__)

    with open(f"{p}/resources/pgkb/action_table.txt") as f:
        next(f)
        for line in f:
            fields = line.strip().split("\t")
            chemical = fields[0]
            gene = fields[1]
            url = fields[2]
            summary = fields[4]
            phenotype = fields[5]
            action = fields[6]

            if chemical not in result:
                result[chemical] = {}

            if gene not in result[chemical]:
                result[chemical][gene] = {
                    "summary": summary,
                    "url": url,
                    "pt": {},
                }

            if phenotype not in result[chemical][gene]["pt"]:
                result[chemical][gene]["pt"][phenotype] = action

    return result

def _read_pair_table():
    result1 = []
    p = os.path.dirname(__file__)

    with open(f"{p}/resources/cpic/cpicPairs.csv") as f:
        result2 = next(f).strip().strip('"').replace(
            "Date last updated: ", "")
        header = next(f).strip().split(",")
        for line in f:
            fields = line.replace(", ", "/").strip().split(",")
            t = [x if x else "-" for x in fields]
            t = dict(zip(header, t))
            result1.append(t)

    return result1, result2

def _read_genotype_file(f, fn, target_genes):
    result1 = {}
    result2 = ""

    if fn:
        f = open(fn)

    header = next(f).strip().split("\t")

    for line in f:
        fields = line.strip().split("\t")
        gene = fields[0]

        if not result2:
            result2 = fields[1]

        result1[gene] = dict(zip(header, fields))

    if fn:
        f.close()

    for target_gene in target_genes:
        if target_gene not in result1:
            result1[target_gene] = dict(zip(header, ["-" for x in header]))

    return result1, result2

def gt2html(
        fn: str,
        f: Optional[TextIO] = None,
        **kwargs
    ) -> str:
    """
    Create HTML report using Stargazer data.

    Returns:
        str: HTML report.

    Args:
        fn (str): Genotype file.
        f (TextIO, optional): Genotype file.
    """

    target_genes = get_target_genes()
    action_table = _read_action_table()
    pair_table, cpic_date = _read_pair_table()
    genotype_table, sample_id = _read_genotype_file(f, fn, target_genes)

    for k, v in genotype_table.items():
        if v["hap1_main"] == "-":
            pt = "-"
        else:
            pt = phenotyper(k, v["hap1_main"], v["hap2_main"])
        genotype_table[k]["phenotype"] = pt

    assessed = [x for x in target_genes if genotype_table[x]["status"] != "-"]
    typed = [x for x in target_genes if genotype_table[x]["status"] == "g"]

    string = (
        "<!DOCTYPE html>\n"
        "<html>\n"
        "<head>\n"
        "<title>Stargazer Report</title>\n"
        "<style>\n"
        "* {\n"
        "  font-family: Arial, Helvetica, sans-serif;\n"
        "}\n"
        "table {\n"
        "  border-collapse: collapse;\n"
        "  width: 100%;\n"
        "  font-size: 80%;\n"
        "}\n"
        "th, td {\n"
        "  border: 1px solid black;\n"
        "  padding: 4px;\n"
        "}\n"
        "</style>\n"
        "</head>\n"
        "<body>\n"
        "<h1>Stargazer Report</h1>\n"
        "<p>\n"
        f"  Sample ID: {sample_id}<br />\n"
        f"  Date: {datetime.datetime.now()}<br />\n"
        f"  Genes examined: {len(assessed)}/{len(target_genes)}<br />\n"
        f"  Genotypes called: {len(typed)}/{len(assessed)}<br />\n"
        "</p>\n"
        "<h2>Introduction</h2>\n"
        "<p>\n"
        "  Thank you for choosing Stargazer! Stargazer is a \n"
        "  bioinformatiscs tool for predicting how a person's DNA \n"
        "  affects their response to hundreds of medications. \n"
        "  Stargazer does this by accurately calling star alleles \n"
        "  (haplotypes) in pharmacogenetic (PGx) genes, which are \n"
        "  defined by single-nucleotide variants (SNVs), small \n"
        "  insertion-deletions (indels), and/or large structural \n"
        "  variants (SVs). Once identified, these star alleles can \n"
        "  be translated to an activity score (AS), which is in turn \n"
        "  used to predict the person's drug response. Stargazer can \n"
        "  utilize genomic data from various sources including \n"
        "  next-generation sequencing (NGS) and single nucleotide \n"
        "  polymorphism (SNP) array. For NGS data, Stargazer supports \n"
        "  whole genome sequencing (WGS) as well as targeted \n"
        "  sequencing such as whole exome sequencing (WES). \n"
        "  For more details please visit the Stargazer website \n"
        "  (https://stargazer.gs.washington.edu/stargazerweb/).\n"
        "</p>\n"
        "<p>\n"
        "  This report includes hundreds of gene/drug pairs \n"
        "  (e.g., CYP2D6/codeine) with accompanying levels of evidence \n"
        "  for changing drug choice and dosing decisions. These pairs \n"
        "  are described by the Clinical Pharmacogenetics Implementation \n"
        "  Consortium (CPIC), which is an international assoication \n"
        "  whose primary goal is to facilitate the use of PGx tests for \n"
        "  patient care. Most importantly, CPIC provides detailed \n"
        "  guidelines for helping clinicians understand how available \n"
        "  genetic test results should be used to optimize drug \n"
        f"  therapy. As of {cpic_date}, there are \n"
        f"  {len(pair_table)} gene/drug pairs listed in the CPIC \n"
        "  website (https://cpicpgx.org). Finally, the Food and Drug \n"
        "  Administration (FDA) provides additional guidance by \n"
        "  requiring applicable PGx test information be included in \n"
        "  the drug labeling. These FDA-approved drug labels are \n"
        "  included in this report as well.\n"
        "</p>\n"
        "<p>\n"
        "  Disclaimer: This report is still very much in development. \n"
        "  Please do not use it other than for code testing. Thank you.\n"
        "</p>\n"
        "<h2>Sections</h2>\n"
        "<ul>\n"
        "  <li>Overview</li>\n"
        "  <li>Genotypes</li>\n"
        "  <li>Drugs</li>\n"
        "  <li>Recommendations</li>\n"
        "</ul>\n"
        "<p style='page-break-before: always;'>\n"
        "<h2>Overview</h2>\n"
        "<p>\n"
        "  PGx genes whose genotype leads to altered phenotype are \n"
        "  shown in <span style='color: red;'>red</span>.\n"
        "</p>\n"
        f"{_add_overview_section(genotype_table, pair_table, target_genes)}\n"
        "<p style='page-break-before: always;'>\n"
        "<h2>Genotypes</h2>\n"
        f"{_add_genotypes_section(genotype_table, target_genes)}\n"
        "<p style='page-break-before: always;'>\n"
        "<h2>Drugs</h2>\n"
        "<p>\n"
        "  Gene/drug pairs are shown in \n"
        "  <span style='font-weight: bold;'>bold</span> if genotype is \n"
        "  available, and in \n"
        "  <span style='font-weight: bold; color: red;'>red</span> \n"
        "  if altered phenotype is predicted.\n"
        "</p>\n"
        f"{_add_drugs_section(genotype_table, pair_table)}\n"
        "<p style='page-break-before: always;'>\n"
        "<h2>Recommendations</h2>\n"
        "<p>\n"
        "  Gene/drug pairs are shown in \n"
        "  <span style='font-weight: bold;'>bold</span> if genotype is \n"
        "  available, and in \n"
        "  <span style='font-weight: bold; color: red;'>red</span> \n"
        "  if altered phenotype is predicted.\n"
        "</p>\n"
        f"{_add_recommendations_section(action_table)}\n"
        "</body>\n"
        "</html>\n"
    )

    return string



        "       <!-- Core theme CSS (includes Bootstrap)-->\n"
        "       <link href='freelancer/css/styles.css' rel='stylesheet'>\n"
        "       <!-- Fonts CSS-->\n"
        "       <link rel='stylesheet' href='freelancer/css/heading.css'>\n"
        "       <link rel='stylesheet' href='freelancer/css/body.css'>\n"


                "<div class='divider-custom-icon'><i class='fas fa-star'></i></div>\n"

                "       <link rel='stylesheet' href='https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css'>\n"



        "       <!-- Font Aweasome -->\n"
        "       <script src='https://kit.fontawesome.com/50a42c0fa1.js' crossorigin='anonymous'></script>\n"

        "<style type='text/css'>\n"
        ".masthead {\n"
        "	padding-top: calc(6rem + 74px);\n"
        "	padding-bottom: 6rem\n"
        "}\n"
        ".masthead .masthead-heading {\n"
        "	font-size: 2.75rem;\n"
        "	line-height: 2.75rem\n"
        "}\n"
        ".masthead .masthead-subheading {\n"
        "	font-size: 1.25rem\n"
        "}\n"
        ".masthead .masthead-avatar {\n"
        "	width: 15rem\n"
        "}\n"
        "@media (min-width:992px) {\n"
        "	.masthead {\n"
        "		padding-top: calc(6rem + 104px);\n"
        "		padding-bottom: 6rem\n"
        "	}\n"
        "	.masthead .masthead-heading {\n"
        "		font-size: 4rem;\n"
        "		line-height: 3.5rem\n"
        "	}\n"
        "	.masthead .masthead-subheading {\n"
        "		font-size: 1.5rem\n"
        "	}\n"
        "}\n"
        ".d-flex {\n"
        "	display: flex!important\n"
        "}\n"
        ".flex-column {\n"
        "	flex-direction: column!important\n"
        "}\n"
        ".mb-5,\n"
        ".my-5 {\n"
        "	margin-bottom: 3rem!important\n"
        "}\n"
        ".mb-0,\n"
        ".my-0 {\n"
        "	margin-bottom: 0!important\n"
        "}\n"
        ".divider-custom {\n"
        "   margin: 1.25rem 0 1.5rem;\n"
        "   width: 100%;\n"
        "   display: flex;\n"
        "   justify-content: center;\n"
        "   align-items: center\n"
        "}\n"
        ".divider-custom .divider-custom-line {\n"
        "   width: 100%;\n"
        "   max-width: 7rem;\n"
        "   height: .25rem;\n"
        "   background-color: #2c3e50;\n"
        "   border-radius: 1rem;\n"
        "   border-color: #2c3e50!important\n"
        "}\n"
        ".divider-custom .divider-custom-line:first-child {\n"
        "   margin-right: 1rem\n"
        "}\n"
        ".divider-custom .divider-custom-line:last-child {\n"
        "   margin-left: 1rem\n"
        "}\n"
        ".divider-custom .divider-custom-icon {\n"
        "   color: #2c3e50!important;\n"
        "   font-size: 2rem\n"
        "}\n"
        ".divider-custom.divider-light .divider-custom-line {\n"
        "   background-color: #fff\n"
        "}\n"
        ".divider-custom.divider-light .divider-custom-icon {\n"
        "   color: #fff!important\n"
        "}\n"
        ".pre-wrap {\n"
        "   white-space: pre-wrap\n"
        "}\n"
        ".font-weight-light {\n"
        "   font-weight: 300!important\n"
        "}\n"
        ".page-section {\n"
	    "padding: 6rem 0\n"
        "}\n"
        ".page-section .page-section-heading {\n"
	    "font-size: 2.25rem;\n"
	    "line-height: 2rem\n"
        "}\n"
        "@media (min-width:992px) {\n"
        "	.page-section .page-section-heading {\n"
        "		font-size: 3rem;\n"
        "		line-height: 2.5rem\n"
        "	}\n"
        "}\n"
        ".page-section-heading {\n"
	    "display: inline-block\n"
        "}\n"
        ".d-inline-block {\n"
	    "display: inline-block!important\n"
        "}\n"
        "@media (min-width:992px) {\n"
        "	.col-lg-4 {\n"
		"            flex: 0 0 33.3333333333%;\n"
		"            max-width: 33.3333333333%\n"
        "   }\n"
        ".col-lg-4 {\n"
	    "position: relative;\n"
	    "width: 100%;\n"
	    "padding-right: 1.5rem;\n"
	    "padding-left: 1.5rem\n"
        "}\n"
        "@media (min-width:992px) {\n"
        "	.col-lg-10 {\n"
		"            flex: 0 0 40%;\n"
		"            max-width: 40%\n"
        "   }\n"
        ".col-lg-10 {\n"
	    "position: relative;\n"
	    "width: 100%;\n"
	    "padding-right: 1.5rem;\n"
	    "padding-left: 1.5rem\n"
        "}\n"
        ".ml-auto,\n"
        ".mx-auto {\n"
	    "        margin-left: auto!important\n"
        "}\n"
        ".lead {\n"
	    "font-size: 1.25rem;\n"
	    "font-weight: 300\n"
        "}\n"
        ".mr-auto,\n"
        ".mx-auto {\n"
	    "margin-right: auto!important\n"
        "}\n"
        ".row {\n"
	    "display: flex;\n"
	    "flex-wrap: wrap;\n"
        "margin-right: -1.5rem;\n"
	    "margin-left: -1.5rem\n"
        "}\n"
        ".justify-content-center {\n"
	    "justify-content: center!important\n"
        "}\n"
        ".align-items-center {\n"
	    "align-items: center!important\n"
        "}\n"
        ".icon-contact {\n"
	    "display: flex;\n"
	    "height: 5rem;\n"
	    "width: 5rem;\n"
        "align-items: center;\n"
        "justify-content: center;\n"
        "font-size: 2rem;\n"
        "background-color: #1abc9c;\n"
        "color: #fff;\n"
        "border-radius: 100%\n"
        "}\n"
        ".text-muted {\n"
        "    color: #6c757d!important\n"
        "}\n"
        ".font-weight-bold {\n"
        "font-weight: 700!important\n"
        "}\n"
        ".font-weight-bolder {\n"
        "    font-weight: bolder!important\n"
        "}\n"
        ".mb-3,\n"
        ".my-3 {\n"
        "    margin-bottom: 1rem!important\n"
        "}\n"
        ".footer {\n"
        "padding-top: 5rem;\n"
        "padding-bottom: 5rem;\n"
        "background-color: #2c3e50;\n"
        "color: #fff\n"
        "}\n"
        "@media (min-width:992px) {\n"
        ".mb-lg-0,\n"
        ".my-lg-0 {\n"
        "    margin-bottom: 0!important\n"
        "}\n"
        ".mb-4,\n"
        ".my-4 {\n"
        "    margin-bottom: 1.5rem!important\n"
        "}\n"
        ".ml-1,\n"
        ".mx-1 {\n"
        "    margin-left: .25rem!important\n"
        "}\n"
        ".mr-1,\n"
        ".mx-1 {\n"
        "    margin-right: .25rem!important\n"
        "}\n"
        ".btn-social {\n"
        "border-radius: 100%;\n"
        "display: inline-flex;\n"
        "width: 3.25rem;\n"
        "height: 3.25rem;\n"
        "font-size: 1.25rem;\n"
        "justify-content: center;\n"
        "align-items: center\n"
        "}\n"
        ".btn-social {\n"
	    "margin-bottom: .5rem\n"
        "}\n"
        ".navbar {\n"
        "    position: relative;\n"
        "    padding: .5rem 1rem\n"
        "}\n"

        ".navbar,\n"
        ".navbar .container,\n"
        ".navbar .container-fluid,\n"
        ".navbar .container-lg,\n"
        ".navbar .container-md,\n"
        ".navbar .container-sm,\n"
        ".navbar .container-xl {\n"
        "    display: flex;\n"
        "    flex-wrap: wrap;\n"
        "    align-items: center;\n"
        "    justify-content: space-between\n"
        "}\n"
        f"{css_string}\n"
        "</style>\n"



        "<script type='text/javascript'>\n"
        "function a.js-scroll-trigger[href*='#']:not([href='#'])').click(function() {\n"
        "    if (location.pathname.replace(/^\//, '') == this.pathname.replace(/^\//, '') && location.hostname == this.hostname) {\n"
        "        var target = $(this.hash);\n"
        "        target = target.length ? target : $('[name=' + this.hash.slice(1) + ']');\n"
        "        if (target.length) {\n"
        "        $('html, body').animate({\n"
        "            scrollTop: (target.offset().top - 71)\n"
        "        }, 1000, 'easeInOutExpo');\n"
        "        return false;\n"
        "        }\n"
        "    }\n"
        "    });\n"
        "</script>\n"
        "<script type='text/javascript'>\n"

        "</script>\n"


        f"      <script src={js_file}></script>\n"

                "       <script src='https://kit.fontawesome.com/50a42c0fa1.js' crossorigin='anonymous'></script>\n"

, shrink-to-fit=no



        "    <nav class='navbar navbar-expand-lg bg-secondary fixed-top' id='mainNav'>\n"
        "        <div class='container'><a class='navbar-brand js-scroll-trigger' href='#page-top'>PGxSnake</a>\n"
        "            <button class='navbar-toggler navbar-toggler-right font-weight-bold bg-primary text-white rounded' type='button' data-toggle='collapse' data-target='#navbarResponsive' aria-controls='navbarResponsive' aria-expanded='false' aria-label='Toggle navigation'>Menu <i class='fas fa-bars'></i></button>\n"
        "            <div class='collapse navbar-collapse' id='navbarResponsive'>\n"
        "                <ul class='navbar-nav ml-auto'>\n"
        "                    <li class='nav-item mx-0 mx-lg-1'><a class='nav-link py-3 px-0 px-lg-3 rounded js-scroll-trigger' href='#Introduction'>INTRODUCTION</a>\n"
        "                    </li>\n"
        "                    <li class='nav-item mx-0 mx-lg-1'><a class='nav-link py-3 px-0 px-lg-3 rounded js-scroll-trigger' href='#Overview'>OVERVIEW</a>\n"
        "                    </li>\n"
        "                    <li class='nav-item mx-0 mx-lg-1'><a class='nav-link py-3 px-0 px-lg-3 rounded js-scroll-trigger' href='#Recommendations'>RECOMMENDATIONS</a>"
        "                    </li>\n"
        "                    <li class='nav-item mx-0 mx-lg-1'><a class='nav-link py-3 px-0 px-lg-3 rounded js-scroll-trigger' href='#Contact'>CONTACT</a>"
        "                    </li>\n"
        "                </ul>\n"
        "            </div>\n"
        "        </div>\n"
        "    </nav>\n"


                "                    <li class='nav-item mx-0 mx-lg-1'><a class='nav-link py-3 px-0 px-lg-3 rounded js-scroll-trigger' href='#Contact'>CONTACT</a>"
        "                    </li>\n"
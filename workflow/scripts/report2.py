import os
import csv
import sys
import datetime
import pkgutil

import pandas as pd
import numpy as np
from pypgx.common import get_target_genes
from pypgx.gt2pt import phenotyper


from recommendations import *

def _read_genotype_file(fn, target_genes):
    result1 = {}
    result2 = ""

    f = open(fn)

    header = next(f).strip().split("\t")

    for line in f:
        fields = line.strip().split("\t")
        gene = fields[0]

        if not result2:
            result2 = fields[1]

        result1[gene] = dict(zip(header, fields))

    
    f.close()

    for target_gene in target_genes:
        if target_gene not in result1:
            result1[target_gene] = dict(zip(header, ["-" for x in header]))

    return result1, result2


def _add_overview_section(genotype_table, target_genes):
    table = (
        "<div class='table-responsive-md'>\n"
        "<table class='table'>\n"
        "   <thead style='text-align:center' class='thead-dark'>"
        "       <tr>\n"
        "           <th scope='col' style='width: 10%;'>No.</th>\n"
        "           <th scope='col' style='width: 30%;'>Gene</th>\n"
        "           <th scope='col' style='width: 30%;'>Genotype</th>\n"
        "           <th scope='col' style='width: 30%;'>Predicted Phenotype</th>\n"
        "       </tr>\n"
        "   </thead>\n"
        "   <tbody>\n"

    )
    normal = ["normal"]
    unknown = ["unknown"]
    nogenotype = ["-"]
    for i, gene in enumerate(target_genes, 1):
        if any([x in genotype_table[gene]["phenotype"] for x in normal]):
            color = "table-success"
        elif any([x in genotype_table[gene]["phenotype"] for x in unknown]):
            color = "table-active"
        elif any([x in genotype_table[gene]["phenotype"] for x in nogenotype]):
            color = "table-dark"
        else:
            color = "table-danger"
        genotype = genotype_table[gene]["hap1_main"] + "/" + genotype_table[gene]["hap2_main"]
        phenotype = genotype_table[gene]["phenotype"]
        table += (
            f"      <tr style='text-align:center' class={color}>\n"
            f"          <th scope='row'>{i}</th>\n"
            f"          <td>{gene.upper()}</td>\n"
            f"          <td>{genotype}</td>\n"
            f"          <td>{phenotype}</td>\n"
            "       </tr>\n"
        )
    return table + "</tbody>\n" + "</table>\n" + "</div>\n"


def _add_recommendations_section(phenotype_table, drug_table):
    table = (
        "<table class='table'>\n"
        "   <thead class='thead-dark'>"
        "       <tr>\n"
        "           <th scope='col' style='width: 20%;'>Drug</th>\n"
        "           <th scope='col' style='width: 80%;'>Recommendation</th>\n"
        "       </tr>\n"
        "   </thead>\n"
        "   <tbody>\n"
    )

    for drug in drug_table:
        recommendation =  eval('{}_recommendation({})'.format(drug, phenotype_table))
        if recommendation == "No recommendation":
            next
        else:
            table += (
                "<tr>\n"
                f"  <td>{drug}</td>\n"
                f"  <td>{recommendation}</td>\n"
                "</tr>\n"
            )
    if table == (
        "<table class='table'>\n"
        "   <thead class='thead-dark'>"
        "       <tr>\n"
        "           <th scope='col' style='width: 20%;'>Drug</th>\n"
        "           <th scope='col' style='width: 80%;'>Recommendation</th>\n"
        "       </tr>\n"
        "   </thead>\n"
        "   <tbody>\n"
    ):
        return(
            "<div style='text-align:center' class='alert alert-success' role='alert'>\n"
            "<h4 class='alert-heading'>Συγχαρητήρια!</h4>\n"
            "<p style='text-align:center'>Σύμφωνα με την φαρμακογονιδιωματική ανάλυση που πραγματοποιήθηκε και τον προβλεπόμενο φαινότυπο για κάθε γονίδιο που εξετάστηκε, δεν υπάρχουν αλλαγές στην δοσολογία φαρμάκων που θα πρέπει να προσέξετε.</p>\n"
            "<hr>\n"
            "<p style='text-align:center' class='mb-0'>Disclaimer: Συμβουλευτείτε τον γιατρό σας.</p>\n"
            "</div>\n"
        )
        #return ("Recommendations are not available for these predicted gene phenotypes.")
    return table + "</tbody>\n" + "</table>"

##### GET SELECTED GENES #####
stargazer_target_genes = get_target_genes()

if snakemake.config["ref"]["build"] == "hg38" and "g6pd" in stargazer_target_genes:
    stargazer_target_genes.remove("g6pd")
if snakemake.config["ref"]["build"] == "hg38" and "gstt1" in stargazer_target_genes:
    stargazer_target_genes.remove("gstt1")    


if snakemake.config["params"]["stargazer"]["target_genes"] == "ALL":
    selected_genes = stargazer_target_genes
else:
    selected_genes = []
    for gene in snakemake.config["params"]["stargazer"]["target_genes"].split(","):
        selected_genes.append(gene.strip().lower())
    for gene in selected_genes:
        if gene not in stargazer_target_genes:
            raise ValueError(f"Unrecognized target gene found: {gene}")

##### DEFINE TARGET GENES #####

target_genes = selected_genes

##### GET PROJECT DIR #####
project_dir=snakemake.config["project_dir"]

##### GET SAMPLES #####
samples = pd.read_table(snakemake.config["samples"]).set_index("sample_name", drop=False)


for sample in samples["sample_name"]:
    GENOTYPE = project_dir + "/" + sample + "/genotypes/" + "/genotypes.txt"
    OUTPUT_REPORT = project_dir + "/" + sample + "/report/" + sample + "_PGxSnake_Report.html"
    OUTPUT_PDF = project_dir + "/" + sample + "/report/" + sample + "_PGxSnake_Report.pdf"
    ##### READ GENOTYPE TABLE  #####
    genotype_table, sample_id = _read_genotype_file(GENOTYPE, target_genes)

    ##### GET ASSESED AND TYPED GENES #####

    assessed = [x for x in target_genes if genotype_table[x]["status"] != "-"]
    typed = [x for x in target_genes if genotype_table[x]["status"] == "g"]


    for k, v in genotype_table.items():
        if v["hap1_main"] == "-":
            pt = "-"
        else:
            pt = phenotyper(k, v["hap1_main"], v["hap2_main"])
        genotype_table[k]["phenotype"] = pt

    phenotype_table = {}
    for gene in typed:
        phenotype_table[gene] = genotype_table[gene]["phenotype"]

    drug_table = ["Azathioprine", "Mercaptopurine", "Thioguanine", "Clopidogrel"]

    css_file = snakemake.config["css_file"]
    with open(css_file, 'r') as f:
        css_string = f.read()
    
    js_file = snakemake.config["js_file"]
    with open(js_file, 'r') as f:
        js_string = f.read()
    
    css_heading = snakemake.config["css_file"]
    with open(css_heading, 'r') as f:
        css_heading_string = f.read()

    css_body = snakemake.config["css_file"]
    with open(css_body, 'r') as f:
        css_body_string = f.read()

    string = (
        "<!doctype html>\n"
        "<html lang='en'>\n"
        "   <head>\n"
        "       <!-- Required meta tags -->\n"
        "       <meta charset='utf-8'>\n"
        "       <meta name='viewport' content='width=device-width, initial-scale=1, shrink-to-fit=no'>\n"
        "       <meta name='description' content=''>\n"
        "       <meta name='author' content=''>\n"
        "       <title>PGxSnake Report</title>\n"
        "       <!-- Font Aweasome -->\n"
        "       <link rel='stylesheet' href='https://use.fontawesome.com/releases/v5.14.0/css/all.css' integrity='sha384-HzLeBuhoNPvSl5KYnjx0BT+WB0QEEqLprO+NBkkk5gbc67FTaL7XIGa2w1L0Xbgc' crossorigin='anonymous'>\n"
        "       <!-- Core theme CSS (includes Bootstrap)-->\n"
        f"      <style type='text/css'>{css_string}</style>\n"
        f"      <style type='text/css'>{css_heading_string}</style>\n"
        f"      <style type='text/css'>{css_body_string}</style>\n"       
        "       <!-- Place this tag in your head or just before your close body tag. -->\n"
        "       <script async defer src='https://buttons.github.io/buttons.js'></script>\n"
        "       <link rel='stylesheet' href='https://unpkg.com/bootstrap-table@1.17.1/dist/bootstrap-table.min.css'>\n"
        "   </head>\n"
        "<body id='page-top'>\n"
        "   <header class='masthead bg-primary text-white text-center'\n"
        "       <div class='container d-flex align-items-center flex-column'>\n"
        "       <!-- Masthead Avatar Image--><img class='masthead-avatar mb-5' src='https://raw.githubusercontent.com/kokyriakidis/Pharmakon/master/bowl-hygieia-symbol-pharmacy-logo_6427-394-removebg-preview.png' alt=''>\n"
        "       <!-- Masthead Heading-->\n"
        "       <h1 class='masthead-heading mb-0'>PGxSnake PGx Report</h1>\n"
        "                <!-- Icon Divider-->\n"
        "       <div class='divider-custom divider-light'>\n"
        "           <div class='divider-custom-line'></div>\n"
        "           <div class='divider-custom-icon'><i class='fas fa-star'></i></div>\n"
        "           <div class='divider-custom-line'></div>\n"
        "       </div>\n"
        "       <!-- Masthead Subheading-->\n"
        "       <p class='pre-wrap masthead-subheading font-weight-light mb-0'>Fast - Accurate - Reproducible</p>\n"
        "       <p>\n</p>\n"
        "       <a href='https://github.com/kokyriakidis?tab=repositories' type='button' class='btn btn-warning' target='_blank' aria-pressed='true' >\n"
        "           <span class='fab fa-github'>\n"
        "           Star on GitHub\n"
        "           </span>\n"
        "       </a>" 
        "       </div>\n"
        "   </header>\n"
        "<section class='footer text-center text-white mb-0'>\n"
        "    <div class='container'>\n"
        "        <!-- Sample Section Heading-->\n"
        "        <div class='text-center'>\n"
        f"            <h2 class='page-section-heading d-inline-block text-white'>Sample ID: {sample}</h2>\n"
        "        </div>\n"
        "        <!-- Icon Divider-->\n"
        "        <div class='divider-custom divider-light'>\n"
        "            <div class='divider-custom-line'></div>\n"
        "            <div class='divider-custom-icon'><i class='fas fa-star'></i></div>\n"
        "            <div class='divider-custom-line'></div>\n"
        "        </div>\n"
        "        <div class='row justify-content-center'>\n"
        "            <!-- Footer Location-->\n"
        "            <div class='col-lg-5 mb-5 mb-lg-0'>\n"
        "                <h4 class='mb-4'>DATE</h4>\n"
        f"                <p class='pre-wrap lead mb-0'>{datetime.datetime.now().date()}</p>\n"
        "            </div>\n"
        "            <!-- Footer Social Icons-->\n"
        "            <div class='col-lg-5 mb-5 mb-lg-0'>\n"
        "                <h4 class='mb-4'>GENES EXAMINED</h4>\n"
        f"                <p class='pre-wrap lead mb-0'>{len(assessed)}/{len(target_genes)}</p>\n"        
        "            </div>\n"
        "            <!-- Footer About Text-->\n"
        "            <div class='col-lg-5'>\n"
        "                <h4 class='mb-4'>GENES GENOTYPED</h4>\n"
        f"                <p class='pre-wrap lead mb-0'>{len(typed)}/{len(assessed)}</p>\n"  
        "            </div>\n"
        "        </div>\n"
        "    </div>\n"
        "</section>\n"
        "<section class='page-section' id='Introduction'>\n"
        "    <div class='container'>\n"
        "        <!-- Introduction Section Heading-->\n"
        "        <div class='text-center'>\n"
        "            <h2 class='page-section-heading text-secondary d-inline-block mb-0'>INTRODUCTION</h2>\n"
        "        </div>\n"
        "        <!-- Icon Divider-->\n"
        "        <div class='divider-custom'>\n"
        "            <div class='divider-custom-line'></div>\n"
        "            <div class='divider-custom-icon'><i class='fas fa-star'></i></div>\n"
        "            <div class='divider-custom-line'></div>\n"
        "        </div>\n"
        "        <!-- Introduction Section Content-->\n"
        "        <div class='row justify-content-center'>\n"
        "            <div class='col-lg-5 ml-auto'>\n"
        "               <p style='text-align:justify' class='pre-wrap lead'>Thank you for choosing PGxSnake! PGxSnake is a bioinformatiscs tool for predicting how a person's DNA affects their response to hundreds of medications. PGxSnake does this by accurately calling star alleles (haplotypes), using Stargazer, in pharmacogenetic (PGx) genes, which are defined by single-nucleotide variants (SNVs), small insertion-deletions (indels), and/or large structural variants (SVs). Once identified, these star alleles can be translated to an activity score (AS), which is in turn used to predict the person's drug response.</p>\n"
        "            </div>\n"
        "            <div class='col-lg-5 mr-auto'>\n"
        "               <p style='text-align:justify' class='pre-wrap lead'>This report includes gene/drug pairs (e.g. CYP2D6/codeine) with accompanying levels of evidence for changing drug choice and dosing decisions. These pairs are described by the Clinical Pharmacogenetics Implementation Consortium (CPIC), which is an international assoication whose primary goal is to facilitate the use of PGx tests for patient care. Most importantly, CPIC provides detailed guidelines for helping clinicians understand how available genetic test results should be used to optimize drug therapy.</p>\n"
        "            </div>\n"
        "        </div>\n"
        "    </div>\n"
        "</section>\n"
        "<section class='page-section' id='Overview'>\n"
        "    <div class='container'>\n"
        "        <!-- Overview Section Heading-->\n"
        "        <div class='text-center'>\n"
        "            <h2 class='page-section-heading text-secondary d-inline-block mb-0'>OVERVIEW</h2>\n"
        "        </div>\n"
        "        <!-- Icon Divider-->\n"
        "        <div class='divider-custom'>\n"
        "            <div class='divider-custom-line'></div>\n"
        "            <div class='divider-custom-icon'><i class='fas fa-star'></i></div>\n"
        "            <div class='divider-custom-line'></div>\n"
        "        </div>\n"
        "        <!-- Introduction Section Content-->\n"
        f"      {_add_overview_section(genotype_table, target_genes)}\n"
        "   </div>\n"
        "</section>\n"
        "<section class='page-section' id='Recommendations'>\n"
        "    <div class='container'>\n"
        "        <!-- Recommendations Section Heading-->\n"
        "        <div class='text-center'>\n"
        "            <h2 class='page-section-heading text-secondary d-inline-block mb-0'>RECOMMENDATIONS</h2>\n"
        "        </div>\n"
        "        <!-- Icon Divider-->\n"
        "        <div class='divider-custom'>\n"
        "            <div class='divider-custom-line'></div>\n"
        "            <div class='divider-custom-icon'><i class='fas fa-star'></i></div>\n"
        "            <div class='divider-custom-line'></div>\n"
        "        </div>\n"
        "        <!-- Introduction Section Content-->\n"
        f"      {_add_recommendations_section(phenotype_table, drug_table)}\n"
        "   </div>\n"
        "</section>\n"
        "<section class='page-section' id='Contact'>\n"
        "    <div class='container'>\n"
        "        <!-- Contact Section Heading-->\n"
        "        <div class='text-center'>\n"
        "            <h2 class='page-section-heading text-secondary d-inline-block mb-0'>CONTACT</h2>\n"
        "        </div>\n"
        "        <!-- Icon Divider-->\n"
        "        <div class='divider-custom'>\n"
        "            <div class='divider-custom-line'></div>\n"
        "            <div class='divider-custom-icon'><i class='fas fa-star'></i></div>\n"
        "            <div class='divider-custom-line'></div>\n"
        "        </div>\n"
        "        <!-- Contact Section Content-->\n"
        "        <div class='row justify-content-center'>\n"
        "            <div class='col-lg-4'>\n"
        "                <div class='d-flex flex-column align-items-center'>\n"
        "                    <div class='icon-contact mb-3'><i class='fas fa-mobile-alt'></i></div>\n"
        "                    <div class='text-muted'>Phone</div>\n"
        "                    <div class='lead font-weight-bold'>(+30) 6981617273</div>\n"
        "                </div>\n"
        "            </div>\n"
        "            <div class='col-lg-4'>\n"
        "                <div class='d-flex flex-column align-items-center'>\n"
        "                    <div class='icon-contact mb-3'><i class='far fa-envelope'></i></div>\n"
        "                    <div class='text-muted'>Email</div><a class='lead font-weight-bold'>kokyriakidis@gmail.com</a>\n"
        "                </div>\n"
        "            </div>\n"
        "        </div>\n"
        "    </div>\n"
        "</section>\n"
        "<footer class='footer text-center'>\n"
        "    <div class='container'>\n"
        "        <div class='row'>\n"
        "            <!-- Footer Location-->\n"
        "            <div class='col-lg-4 mb-5 mb-lg-0'>\n"
        "                <h4 class='mb-4'>LOCATION</h4>\n"
        "                <p class='pre-wrap lead mb-0'>Lab. of Pharmacology, School of Pharmacy, Aristotle University of Thessaloniki, Greece</p>\n"
        "            </div>\n"
        "            <!-- Footer Social Icons-->\n"
        "            <div class='col-lg-4 mb-5 mb-lg-0'>\n"
        "                <h4 class='mb-4'>AROUND THE WEB</h4><a class='btn btn-outline-light btn-social mx-1' href='https://www.facebook.com/kokyriakidis'><i class='fab fa-fw fa-facebook-f'></i></a><a class='btn btn-outline-light btn-social mx-1' href='https://www.twitter.com/kokyriakidis'><i class='fab fa-fw fa-twitter'></i></a><a class='btn btn-outline-light btn-social mx-1' href='https://www.linkedin.com/in/konstantinos-kyriakidis-3b393672'><i class='fab fa-fw fa-linkedin-in'></i></a>\n"
        "            </div>\n"
        "            <!-- Footer About Text-->\n"
        "            <div class='col-lg-4'>\n"
        "                <h4 class='mb-4'>ABOUT PGxSnake</h4>\n"
        "                <p class='pre-wrap lead mb-0'>PGxSnake is a free to use, GPLv3 licensed Pharmacogenomics pipeline created by Konstantinos Kyriakidis</p>\n"
        "            </div>\n"
        "        </div>\n"
        "    </div>\n"
        "</footer>\n"
        "       <!-- Bootstrap core JS-->\n"
        "       <script src='https://cdnjs.cloudflare.com/ajax/libs/jquery/3.4.1/jquery.min.js'></script>\n"
        "       <script src='https://stackpath.bootstrapcdn.com/bootstrap/4.4.1/js/bootstrap.bundle.min.js'></script>\n"
        "       <script src='https://unpkg.com/bootstrap-table@1.17.1/dist/bootstrap-table.min.js'></script>\n"
        "       <!-- Third party plugin JS-->\n"
        "       <script src='https://cdnjs.cloudflare.com/ajax/libs/jquery-easing/1.4.1/jquery.easing.js'></script>\n"
        "       <!-- Core theme JS-->\n"
        f"      <script type='text/javascript'>{js_string}</script>\n"
        "   </body>\n"
        "</html>\n"
    )

    with open(
        OUTPUT_REPORT, "w"
    ) as f:
        f.write(string)































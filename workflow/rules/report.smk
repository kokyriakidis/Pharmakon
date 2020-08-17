
rule pgx_report:
    input:
        expand("{project_dir}/{sample}/genotypes/genotypes.txt", project_dir=config["project_dir"], sample=samples.index)
    output:
        expand("{project_dir}/{sample}/report/{sample}_PGxSnake_Report.html", project_dir=config["project_dir"], sample=samples.index)
    conda:
        "../envs/report.yaml"
    script:
        "../scripts/report2.py"
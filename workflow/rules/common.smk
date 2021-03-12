import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

min_version("5.18.0")

report: "../report/workflow.rst"

container: "continuumio/miniconda3:4.8.2"

###### Config file #####
configfile: "config/config.yaml"

#validate(config, schema="../schemas/config.schema.yaml")

##### Sample sheet #####
samples = pd.read_table(config["samples"]).set_index("samples", drop=False)
#validate(samples, schema="../schemas/samples.schema.yaml")

##### Wildcard constraints #####
wildcard_constraints:
    vartype="snvs|indels",
    sample="|".join(samples.index)


##### Helper functions #####

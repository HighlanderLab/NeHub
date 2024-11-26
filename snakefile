configfile: "./config/global_config.yaml"
import glob
import os
import json
# Define input and output file patterns
workDir=config["workdir"]
vcf=config["vcffile"]
_, vcf_f = os.path.split(vcf)

include: "rules/sampling_vcf.smk"
for method in config["Method"]:
    include: f"rules/r{method}.smk"

# MSMC's covers up to two groups
# SMCpp covers up to two groups
# chromosome must be named chr* in stead of number, otherwise it will cause the failure of SMCpp
rule all:
    input:
        expand(f'{workDir}/{vcf_f}/{{Method}}_results.json', Method=config["Method"])




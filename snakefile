configfile: "./config/global_config.yaml"
import glob
import os
import json
# Define input and output file patterns
workDir=config["workdir"]
vcf=config["vcffile"]
_, vcf_f = os.path.split(vcf)

if "GONE" in config["Method"]:
    include: "rules/sampling_vcf.smk",
    include: "rules/rGONE.smk"

if 'PSMC' in config["Method"]:
    include: "rules/rPSMC.smk"

rule all:
    input:
        expand(f'{workDir}/{vcf_f}/{{Method}}_results.json', Method=config["Method"])



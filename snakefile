# Import the config.yaml file
configfile: "./config/GONE_config.yaml"

# Define input and output file patterns
workDir=config["workdir"]
data=config["vcffile"]
sample_name=config["sampleList"]
plink_path=config["plinkModule"]

# Rule to define all outputs for the workflow
rule all:
    input:
        expand(f'{workDir}/sampling{{number}}.txt', number=range(config["numSampling"]))

rule sampling:
    input:
        sample_name
    output:
        expand(f'{workDir}/sampling{{number}}.txt', number=range(config["numSampling"]))

    conda: "HLab_tsinfer"
    params:
        sampling_num=config["numSampling"],
        sample_size=config["sampleSize"]
    shell:  
        '''
        python3 ./script/sampling.py -i {input} -n {params.sampling_num} -s {params.sample_size} -o {workDir}
        '''

# rule all:
#     input:
#         expand("{data}/{sample_name}.bed", data=data, sample_name=sample_name),
#         expand("{data}/{sample_name}.bim", data=data, sample_name=sample_name),
#         expand("{data}/{sample_name}.fam", data=data, sample_name=sample_name)

# rule run_gone:
#     input:
#         "GONE_config.yaml"
#     output:
#         "GONE_config.yaml"
#     shell:
#         "python3 {configfile['workdir']}/GONE.py {configfile['workdir']}/GONE_config.yaml"
        

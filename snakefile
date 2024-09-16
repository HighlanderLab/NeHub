# Import the config.yaml file
configfile: "GONE_config.yaml"

# Define input and output file patterns
data=configfile["vcffile"]
sample_name=configfile["sampleList"]
plink_path=configfile["plinkModule"]
workdir=configfile["workdir"]

rule sampling:
    input:
        sample_name
    output:
        workdir
    conda: "HLab_tsinfer"
    params:
        sampling_num=configfile["numSampling"],
        sample_size=configfile["sampleSize"]
    shell:  
        '''
        python3 ./script/sampling.py -i {input} -n {params.sampling_num} -s {params.sample_size} -o {output}
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
        

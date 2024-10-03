# Import the config.yaml file
configfile: "./config/GONE_config.yaml"
import glob
import os
# Define input and output file patterns
workDir=config["workdir"]
vcf=config["vcffile"]
_, vcf_f = os.path.split(vcf)


rule all:
    input:
        expand(f'{workDir}/{vcf_f}/{{Method}}_results/sampling{{num}}/Output_Ne_data', num=range(config["numSampling"]),Method=config["Method"])


rule sampling:
    input:
        f'{vcf}.vcf.gz',
        f'{vcf}.popLabel'
    output:
        expand(f'{workDir}/{vcf_f}/sampling{{number}}.vcf.gz', number=range(config["numSampling"]))
    conda: "HLab_tsinfer"
    params:
        bcftools=config["bcftoolsModule"],
        sampling_num=config["numSampling"],
        sample_size=config["sampleSize"],
        outputDir=lambda wildcards: f"{workDir}/{vcf_f}"
    shell:  
        '''
        module load {params.bcftools}
        if [ ! -d {params.outputDir} ]; then
        mkdir -p "{params.outputDir}"
        fi

        python3 ./script/sampling.py -i {input[1]} -n {params.sampling_num} -s {params.sample_size} -o {params.outputDir}

        # Find all sample*.txt files and run bcftools view on each
        for sample_file in {params.outputDir}/sampling*.txt; do
            base_name=$(basename "$sample_file" .txt)
            output_vcf="{params.outputDir}/${{base_name}}.vcf.gz"
            # Run bcftools view to extract the sample from the VCF
            bcftools view -S "${{sample_file}}" {input[0]} -o $output_vcf
        done
        '''

rule run_gone:
    input:
        f'{workDir}/{vcf_f}/{{file_number}}.vcf.gz'
    output:
        # Expanding OUTPUT_data for each sampling number
        expand(f'{workDir}/{vcf_f}/GONE_results/{{file_number}}/Output_Ne_data',file_number='{file_number}')
    params:
        plink_path=config["plinkModule"],
        GONE=config["GONE_path"],  # GONE path from config
        GONE_setting=config["GONE_setting"]
        work_dir=lambda wildcards: f'{workDir}/{vcf_f}/'
    shell:
        '''
        module load {params.plink_path}
        cd {params.work_dir}
        mkdir -p GONE_results
        cd GONE_results
        parallel -j {threads} '
            base_name=$(basename {{}} .vcf.gz)
            mkdir -p "$base_name"
            cd "$base_name"
            plink --vcf {{}} --recode --out data
            # Create symbolic links
            ln -sf {params.GONE}Linux/script_GONE.sh ./script_GONE.sh
            ln -sf {params.GONE}Linux/PROGRAMMES/ ./
            ln -sf {params.GONE_setting} ./INPUT_PARAMETERS_FILE
            # Run the GONE script
            time bash script_GONE.sh data > "GONE.stdout" 2> "GONE.stderr"
        ' ::: {input}
        '''

# rule summary_GONE:
#     input:
#     output:
#     params:
#     shell:
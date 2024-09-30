# Import the config.yaml file
configfile: "./config/GONE_config.yaml"

# Define input and output file patterns
workDir=config["workdir"]
sample_name=config["sampleList"]
vcffile=config["vcffile"]
vcfdir=config["vcfdir"]
new_work_folder=vcffile.split(".")[0]

rule all:
    input:
        expand(f'{workDir}/{new_work_folder}/sampling{{num}}/OUTPUT_data', num=range(config["numSampling"]))


rule sampling:
    input:
        sample_name
    output:
        expand(f'{workDir}/{new_work_folder}/sampling{{number}}.txt', number=range(config["numSampling"]))
    conda: "HLab_tsinfer"
    params:
        sampling_num=config["numSampling"],
        sample_size=config["sampleSize"],
        outputDir=lambda wildcards: f"{workDir}/{new_work_folder}/"
    shell:  
        '''
        if [ ! -d {params.outputDir} ]; then
        mkdir -p "{params.outputDir}"
        fi
        python3 ./script/sampling.py -i {input} -n {params.sampling_num} -s {params.sample_size} -o {params.outputDir}
        '''

rule prepare_input:
    input:
        rules.sampling.output
    output:
        expand(f'{workDir}/{new_work_folder}/sampling{{num}}/data.log', num=range(config["numSampling"]))
    params:
        plink_path=config["plinkModule"],
        vcf=lambda wildcards: f"{vcfdir}/{vcffile}",
        outputDir=lambda wildcards: f"{workDir}/{new_work_folder}/"
    shell:
        '''
        module load {params.plink_path}
        cd {params.outputDir}

        # Loop over each input file and process it
        for file in {input}; do
            folder_name=$(basename "$file" .txt)  # Remove '.txt' to create folder name
            mkdir -p "$folder_name"  # Create directory based on input file
            plink --vcf {params.vcf} --keep "$file" --recode --out "$folder_name/data"
        done
        '''

rule run_gone:
    input:
        rules.prepare_input.output  # Reference the output of prepare_input rule
    output:
        # Expanding OUTPUT_data for each sampling number
        expand(f'{workDir}/{new_work_folder}/sampling{{num}}/OUTPUT_data', num=range(config["numSampling"]))
    params:
        GONE=config["GONE_path"],  # GONE path from config
        vcfdir=lambda wildcards: vcfdir
    shell:
        '''
        parallel -j {threads} '
            folder_name=$(dirname {{}})  # Get the folder for each file
            cd "$folder_name"

            # Create symbolic links
            ln -sf {params.vcfdir}/script_GONE.sh ./script_GONE.sh
            ln -sf {params.GONE}Linux2/PROGRAMMES/ ./
            ln -sf {params.vcfdir}/INPUT_PARAMETERS_FILE ./INPUT_PARAMETERS_FILE

            # Run the GONE script
            time bash script_GONE.sh data > "GONE.stdout" 2> "GONE.stderr"
        ' ::: {input}
        '''
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
    shell:  
        '''
        if [ ! -d {workDir}/{new_work_folder} ]; then
        mkdir -p "{workDir}/{new_work_folder}"
        fi
        python3 ./script/sampling.py -i {input} -n {params.sampling_num} -s {params.sample_size} -o {workDir}/{new_work_folder}/
        '''

rule prepare_input:
    input:
        rules.sampling.output
    output:
        expand(f'{workDir}/{new_work_folder}/sampling{{num}}/data.log', num=range(config["numSampling"]))
    params:
        plink_path=config["plinkModule"]
    shell:
        '''
        module load {params.plink_path}
        cd {workDir}/{new_work_folder}

        # Loop over each input file and process it
        for file in {input}; do
            folder_name=$(basename "$file" .txt)  # Remove '.txt' to create folder name
            mkdir -p "$folder_name"  # Create directory based on input file
            plink --vcf {vcfdir}/{vcffile} --keep "$file" --recode --out "$folder_name/data"
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
    log: 
        # Log file path
        f'{workDir}/{new_work_folder}/GONE.log'
    shell:
        '''
        for file in {input}; do
            folder_name=$(dirname "$file")  # Get the folder for each file
            cd "$folder_name"

            # Create symbolic links
            ln -sf {vcfdir}/script_GONE.sh ./script_GONE.sh
            ln -sf {params.GONE}Linux2/PROGRAMMES/ ./
            ln -sf {vcfdir}/INPUT_PARAMETERS_FILE ./INPUT_PARAMETERS_FILE

            # Run the GONE script
            time bash script_GONE.sh data
        done
        '''
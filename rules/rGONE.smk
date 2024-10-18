configfile: "./config/GONE_config.yaml"

rule run_gone:
    input:
        f'{workDir}/{vcf_f}/{{file_number}}.vcf.gz'
    output:
        # Expanding OUTPUT_data for each sampling number
        expand(f'{workDir}/{vcf_f}/GONE_results/{{file_number}}/Output_Ne_data',file_number='{file_number}')
    params:
        plink_path=config["plinkModule"],
        GONE=config["GONE_path"],  # GONE path from config
        GONE_setting=config["GONE_setting"],
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

rule summary_GONE:
    input:
        expand(f'{workDir}/{vcf_f}/GONE_results/sampling{{num}}/Output_Ne_data', num=range(config["numSampling"]))
    output:
        f'{workDir}/{vcf_f}/GONE_results.json'
    run:
        combined_li = {}
        for input_file in input:
            with open(str(input_file),"r") as f:
                annotation= f.readline()
                header= f.readline()
                line= f.readline()
                dic={}
                while line !="":
                    li = line.strip("\n").split('\t')
                    dic[int(li[0])]=float(li[1])
                    line= f.readline()
                combined_li[input_file]=dic
        with open(output[0], 'w') as json_file:
            json.dump(combined_li, json_file, indent=4)
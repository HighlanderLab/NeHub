configfile: "./config/MSMC_config.yaml"
import re 
import glob
import os
rule sampling_MSMC:
    input:
        f'{vcf}.vcf.gz',
    output:
        f'{workDir}/{vcf_f}/msmc_results/{{sample}}/finish_flag'
    params:
        popLabel = f'{vcf}.popLabel',
        workDir = f"{workDir}/{vcf_f}/msmc_results",
        num_sampling = config["num_sampling_msmc"],
        sample_size = config["sample_size_msmc"]
    conda: "HLab_tsinfer"
    shell:
        '''
        if [ ! -d {params.workDir} ]; then
        mkdir -p "{params.workDir}"
        fi
        python ./script/sampling_msmc.py  -i {input[0]} -p {params.popLabel} -w {params.workDir} -n {params.num_sampling} -s {params.sample_size}  
        '''

rule run_MSMC:
    input: 
        f'{workDir}/{vcf_f}/msmc_results/{{sample}}/finish_flag'
    output:
        f'{workDir}/{vcf_f}/msmc_results/{{sample}}/msmc_0.final.txt'
    conda: "ne_pipe_msmc2"
    threads: 4
    params:
        tsp = config["timeSegmentPattern"],
        roM=config["rhoOverMu"],
        it=config["maxIterative"]
    shell:
        '''
        parallel -j {threads} '
        dir_name=$(dirname {{}})
        cd $dir_name
        sampling_files=$(ls msmc_*.tab | tr "\\n" " ")
        # Check if the file "finish_flag" exists and read its lines into an array
        count=0
        if [ -f "./finish_flag" ]; then
            while IFS= read -r row; do               
            p=$(echo "$row" | tr -d "\\n")
            echo "$p"
            echo "$sampling_files"
            msmc2_Linux -I $p -i {params.it} -r {params.roM} -p {params.tsp}  -o "msmc_$count" $sampling_files
            count=$((count + 1))
            done < ./finish_flag
        fi
        ' ::: {input}
        '''
        
rule summary_MSMC:
    input:
        expand(f'{workDir}/{vcf_f}/msmc_results/sampling{{sample}}/msmc_0.final.txt', sample=range(config["num_sampling_msmc"]))
    output:
        f'{workDir}/{vcf_f}/MSMC_results.json'
    params:
        mutation_rate = float(config["mutation_rate"])
    run:
        combined_li={}
        for file_path in input:
            directory_name = os.path.dirname(file_path)
            sampling_files = glob.glob(os.path.join(directory_name,"*.final.txt"))
            print(sampling_files)
            li_file ={}
            for i in sampling_files:
                with open(i,"r") as f:
                    content = f.readlines()
                Ne=[]
                for j in content[1:]:
                    li_j=j.strip('\n').split('\t')
                    t=float(li_j[1])
                    la=float(li_j[3])
                    Ne.append([round(t / params.mutation_rate),(1.0 / la) / (2.0 * params.mutation_rate)])
                li_file[os.path.basename(i)]=dict(Ne)

            combined_li[file_path] = li_file

        with open(output[0], 'w') as json_file:
            json.dump(combined_li, json_file, indent=4)

 
 




        




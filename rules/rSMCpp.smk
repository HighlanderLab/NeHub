configfile: "./config/SMCpp_config.yaml"
import re 
import glob
import os
import json

rule prepare_SMCpp:
    input:
        f'{workDir}/{vcf_f}/{{sampling}}.popLabel'
    output:
        f'{workDir}/{vcf_f}/SMCpp_results/{{sampling}}/params_pop.txt'
    run:
        # Inputs
        poplabel_file = input  # Replace with your VCF file
        for i in poplabel_file:
            base_directory = os.path.dirname(i)
            base_name = os.path.basename(i)
            base_name = base_name.replace(".popLabel","")
            # Prepare directories
            os.chdir(base_directory)
            if os.path.isdir("./SMCpp_results"):
                pass
            else:
                os.mkdir("./SMCpp_results")

            os.chdir("./SMCpp_results")
            if os.path.isdir(base_name):
                pass
            else:
                os.mkdir("./{}".format(base_name))
            os.chdir("./{}".format(base_name))

            # Ensure the popLabel file exists
            if not os.path.isfile(i):
                raise FileNotFoundError(f"{i} not found!")

            dic_poplabel={}
            with open(i,"r") as f:
                title=f.readline()
                line = f.readline()
                while line != "":
                    li=line.strip("\n").split("\t")
                    if li[1] not in dic_poplabel:
                        dic_poplabel[li[1]]=[li[0]]
                    else:
                        dic_poplabel[li[1]].append(li[0])
                    line = f.readline()
            with open("params_pop.txt","w") as f1:
                for j in dic_poplabel:
                    f1.write("pop{}:{}\n".format(j,",".join(dic_poplabel[j])))

rule run_SMCpp:
    input:
        f'{workDir}/{vcf_f}/{{sampling}}.vcf.gz',
        f'{workDir}/{vcf_f}/SMCpp_results/{{sampling}}/params_pop.txt'
    output:
        f'{workDir}/{vcf_f}/SMCpp_results/{{sampling}}/results.csv'
    params:
        contig = config["contig"],
        mu=config["mutation_rate"],
        workDir=config["workdir"]
    conda: "bcftools"
    retries: 3
    shell:
        '''
        source /etc/profile.d/modules.sh 
        module load singularity
        parallel -j {threads} '
        # prepare the data and pop samples name
        base_name=$(basename {{}} .vcf.gz)
        directory_name=$(dirname {{}})
        cd $directory_name/SMCpp_results/$base_name
        tabix -f {{}}
        # run SMCpp
        line_count=$(wc -l < params_pop.txt)
        mapfile -t lines < params_pop.txt
        mkdir data_{params.contig}
        if [ "$line_count" -eq 1 ]; then
        singularity run --bind {params.workDir} docker://terhorst/smcpp:latest vcf2smc {{}} data_{params.contig}/pop1.smc.gz {params.contig} ${{lines[0]}}
        singularity run --bind {params.workDir} docker://terhorst/smcpp:latest estimate -o pop1_{params.contig}/ {params.mu} data_{params.contig}/pop1.smc.gz
        singularity run --bind {params.workDir} docker://terhorst/smcpp:latest plot -c results pop1_{params.contig}/model.final.json 

        echo "successfully finish" > finined_flag
        else

        singularity run --bind {params.workDir} docker://terhorst/smcpp:latest vcf2smc {{}} data_{params.contig}/pop1.smc.gz {params.contig} ${{lines[0]}}
        singularity run --bind {params.workDir} docker://terhorst/smcpp:latest estimate -o pop1_{params.contig}/ {params.mu} data_{params.contig}/pop1.smc.gz
        singularity run --bind {params.workDir} docker://terhorst/smcpp:latest vcf2smc {{}} data_{params.contig}/pop2.smc.gz {params.contig} ${{lines[1]}}
        singularity run --bind {params.workDir} docker://terhorst/smcpp:latest estimate -o pop2_{params.contig}/ {params.mu} data_{params.contig}/pop2.smc.gz

        first_row=$(IFS=" "; echo "${{lines[*]}}")
        second_row=$(IFS=" and "; echo "${{lines[*]}}" | tac -s$'\\n')
        singularity run --bind {params.workDir} docker://terhorst/smcpp:latest vcf2smc {{}} data_{params.contig}/pop12.smc.gz {params.contig} $first_row
        singularity run --bind {params.workDir} docker://terhorst/smcpp:latest vcf2smc {{}} data_{params.contig}/pop21.smc.gz {params.contig} $second_row
        singularity run --bind {params.workDir} docker://terhorst/smcpp:latest split -o split_{params.contig}/ pop1_{params.contig}/model.final.json pop2_{params.contig}/model.final.json data_{params.contig}/*.smc.gz
        
        singularity run --bind {params.workDir} docker://terhorst/smcpp:latest plot -c results ./*/model.final.json 

        fi
        ' ::: {input[0]}
        '''

        
        
        
rule summary_SMCpp:
    input:
        expand(f'{workDir}/{vcf_f}/SMCpp_results/sampling{{sample}}/results.csv', sample=range(config["numSampling"]))
    output:
        f'{workDir}/{vcf_f}/SMCpp_results.json'
    run:
        # Define the directory to search
        search_directory = input
        # List to store found files
        combined_li = {}

        # Walk through the directory and subdirectories
        for i in search_directory:
            f=open(i,"r")
            title=f.readline()
            line=f.readline()
            result={}
            while line != "":
                li=line.split(",")
                if li[0] not in result:
                    result[li[0]]={li[1]:li[2]}
                else:
                    result[li[0]][li[1]] = li[2]
                line=f.readline()

            combined_li[i]=result

        with open(output[0], 'w') as json_file:
            json.dump(combined_li, json_file, indent=4)

 
 




        




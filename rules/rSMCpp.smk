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
        vcf=f'{workDir}/{vcf_f}/{{sampling}}.vcf.gz',
        param=f'{workDir}/{vcf_f}/SMCpp_results/{{sampling}}/params_pop.txt'
    output:
        f'{workDir}/{vcf_f}/SMCpp_results/{{sampling}}/pop1/model.final.json'
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
        # prepare the data and pop samples name
        base_name=$(basename {input.vcf} .vcf.gz)
        directory_name=$(dirname {input.vcf})
        cd $directory_name/SMCpp_results/$base_name
        bcftools index -c -f {input.vcf}
        # run SMCpp
        line_count=$(wc -l < params_pop.txt)
        mapfile -t lines < params_pop.txt
        mkdir -p data_{params.contig}
        singularity run --bind {params.workDir} docker://terhorst/smcpp:latest vcf2smc {input.vcf} data_{params.contig}/pop1.smc.gz {params.contig} ${{lines[0]}}
        singularity run --bind {params.workDir} docker://terhorst/smcpp:latest estimate -v -o pop1/ {params.mu} data_{params.contig}/pop1.smc.gz
        '''


rule run_SMCpp_plot:
    input:
        f'{workDir}/{vcf_f}/SMCpp_results/{{sampling}}/pop1/model.final.json'
    output:
        f'{workDir}/{vcf_f}/SMCpp_results/{{sampling}}/results.csv'
    params:
        contig = config["contig"],
        workDir=config["workdir"]
    conda: "bcftools"
    retries: 3
    shell:
        '''
        directory_name=$(dirname {input[0]})
        cd $directory_name
        cd ../
        source /etc/profile.d/modules.sh 
        module load singularity
        singularity run --bind {params.workDir} docker://terhorst/smcpp:latest plot -c results {input[0]}
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

 
 




        




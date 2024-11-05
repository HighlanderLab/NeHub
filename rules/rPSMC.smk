configfile: "./config/PSMC_config.yaml"
import re 

def get_samples(vcf):
    with open(f'{vcf}.popLabel') as f:
        lines = f.readlines()[1:]  # Skip the header
        samples = [line.split("\t")[0] for line in lines]
    return samples
    
rule extract_sample:
    input:
        f'{vcf}.vcf.gz',
        f'{vcf}.popLabel',
        f'{config["ref_genome"]}'
    output:
        expand(f'{workDir}/{vcf_f}/psmc_results/{{sample}}.fq.gz', sample=get_samples(vcf))
    conda: "bcftools"
    shell:
        '''
        mkdir -p {workDir}/{vcf_f}/psmc_results/
        tabix -f {input[0]}
        samples=$(awk 'NR >1 {{print $1}}' {input[1]})
        for sample in ${{samples}}; do
        cat {input[2]} |vcf-consensus -i {input[0]} -s ${{sample}}|\
            sed 's/>/@/g' | bgzip > {workDir}/{vcf_f}/psmc_results/${{sample}}.fq.gz
        done
        '''

rule run_psmc:
    input: 
        f'{workDir}/{vcf_f}/psmc_results/{{sample}}.fq.gz'
    output:
        f'{workDir}/{vcf_f}/psmc_results/{{sample}}.psmc'
    conda: "ne_pipe_psmc"
    params:
        work_dir=lambda wildcards: f'{workDir}/{vcf_f}/psmc_results',
        p = config["psmc_setting"]
    shell:
        '''
        cd {params.work_dir}
        parallel -j {threads} '
        base_name=$(basename {{}} .fq.gz)
        fq2psmcfa -q0 {{}} > ${{base_namev}}.psmcfa
        psmc -N25 -t20 -r5 -p {params.p} ${{base_namev}}.psmcfa -o ${{base_name}}.psmc
        ' ::: {input}
        '''

rule summary_PSMC:
    input:
        expand(f'{workDir}/{vcf_f}/psmc_results/{{sample}}.psmc', sample=get_samples(vcf))
        
    output:
        f'{workDir}/{vcf_f}/PSMC_results.json'
    params:
        bin_size = float(config["bin_size"]),
        mutation_rate = float(config["mutation_rate"])
    run:
        combined_li={}
        for file in input:
            with open(file,"r") as f:
                content = f.read()

                # Split the content by the '//' delimiter
                sections = content.split('//\n')[-2]
            theta_pattern = r"^TR\s+([\d.]+)\s+([\d.]+)"
            rs_pattern = r"^RS\s+\d+\s+([\d.]+)+\s+([\d.]+)"
            Ne=[]
            for i in sections.split("\n"):
                match= re.match( theta_pattern ,i)
                if match:
                    theta = float(match.group(1))
                    print({params.mutation_rate},{params.bin_size})
                    N_0=theta/(4*params.mutation_rate)/params.bin_size
                match=re.match(rs_pattern,i)
                if match:
                    t_k = float(match.group(1))
                    lambda_k = float(match.group(2))
                    try:
                        T_k = 2*N_0 * t_k
                        N_k = N_0 * lambda_k
                    except:
                        print("The .psmc file has issue, please make sure the order of TR line and RS line ")
                    Ne.append([T_k,N_k])            
            combined_li[file]=dict(Ne)
        with open(output[0], 'w') as json_file:
            json.dump(combined_li, json_file, indent=4)






        




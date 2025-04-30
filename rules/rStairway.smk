configfile: "./config/Stairway_config.yaml"

rule polarize:
    input:
        vcf = f'{workDir}/{vcf_f}/sampling{{num}}.vcf.gz',
    output:
        f'{workDir}/{vcf_f}/Stairway_results/sampling{{num}}_polarized.vcf.gz'
    conda:
        "HLab_tsinfer"
    params:
        ancestral = config["ancestral_genome"]
    shell:
        '''
        python ./script/polarize.py -v {input.vcf} -a {params.ancestral} -o {output}
        '''

rule filter_missing:
    input:
        vcf = f'{workDir}/{vcf_f}/Stairway_results/sampling{{num}}_polarized.vcf.gz',
    output:
        f'{workDir}/{vcf_f}/Stairway_results/sampling{{num}}_polarized_filter.vcf.gz'
    conda:
        "bcftools"
    shell:
        '''
        bcftools view -e 'FILTER="MISSING"' {input.vcf} -Oz -o {output}
        '''

rule generate_sfs:
    input:
        vcf = f'{workDir}/{vcf_f}/Stairway_results/sampling{{num}}_polarized_filter.vcf.gz',
        popfile = f'{workDir}/{vcf_f}/sampling{{num}}.popLabel',
    output:
        f'{workDir}/{vcf_f}/Stairway_results/sampling{{num}}_SFS.obs'
    conda:
        "HLab_tsinfer"
    params:
        easySFS=config["easysfs_path"],
        ploid=config["ploids"],
        fold=config["fold"],
    shell:
        '''
        directory_name=$(dirname {input.popfile})
        base_name=$(basename {input.popfile} .popLabel)
        cd $directory_name
        mkdir -p Stairway_results/
        cd Stairway_results/
        # get the the minimum sample number in populations
        awk 'NR > 1' {input.popfile}  > ${{base_name}}.poplabels
        min_count=$(cut -d $'\t' -f2 ${{base_name}}.poplabels | sort | uniq -c | sort -n | head -n 1 | awk '{{print $1}}')
        
        if [ {params.fold} -eq 1 ]; then
            project=$((min_count))
            {params.easySFS} -i {input.vcf} -p ${{base_name}}.poplabels -a --proj $project -o ${{base_name}}_SFS
            for file in ${{base_name}}_SFS/fastsimcoal2/*_MAFpop0.obs; do
                base=$(basename "$file" _MAFpop0.obs)
                new_file="${{base_name}}_SFS.obs"
                cp "$file" "$new_file"
            done
        else
            project=$(({params.ploid} * min_count))
            {params.easySFS} -i {input.vcf} -p ${{base_name}}.poplabels -a --unfold --proj $project -o ${{base_name}}_SFS
            for file in ${{base_name}}_SFS/fastsimcoal2/*_DAFpop0.obs; do
                base=$(basename "$file" _DAFpop0.obs)
                new_file="${{base_name}}_SFS.obs"
                cp "$file" "$new_file"
            done
        fi
        rm -rf ${{base_name}}_SFS
        '''

rule write_stairway_blueprint:
    input:
        vcf = f'{workDir}/{vcf_f}/sampling{{num}}.vcf.gz',
        SFS = f'{workDir}/{vcf_f}/Stairway_results/sampling{{num}}_SFS.obs',
    output:
        blueprint = f'{workDir}/{vcf_f}/Stairway_results/sampling{{num}}.blueprint',
    conda:
        "bcftools"
    params:
        fold=config["fold"],
        mu=config["mutation_rate"],
        ploids=config["ploids"],
        generation_interval=config["generation_interval"],
        stairway_path=config["stairway_config"]
    shell:
        '''
        script_path=$(realpath ./script/write_stairway_blueprint.py)
        directory_name=$(dirname {input.SFS})
        cd $directory_name
        nind=$(bcftools query -l {input.vcf} | wc -l)
        chr=$(bcftools view -h {input.vcf} | grep "^##contig" | awk -F'[,=]' '{{print $3}}')
        Length=$(bcftools view -h {input.vcf} | grep "^##contig" | awk -F'[,=]' '{{print $5}}' | sed 's/>//')

        python $script_path \
            --input {input.SFS} --nind $nind \
            --stairway_plot_dir {params.stairway_path} \
            --chr $chr --length $Length \
            --mutation_rate {params.mu} --fold {params.fold} \
            --ploid {params.ploids} --generation_time {params.generation_interval} \
            --output {output.blueprint}
        '''


rule run_stairway:
    input: 
        f'{workDir}/{vcf_f}/Stairway_results/sampling{{num}}.blueprint',
    output:
        f'{workDir}/{vcf_f}/Stairway_results/sampling{{num}}/sampling{{num}}.final.summary',
    params:
        stairway_path=config["stairway_config"]
    shell:
        '''
        directory_name=$(dirname {input})
        base_name=$(basename {input})
        cd $directory_name
        java -cp {params.stairway_path} Stairbuilder {input}
        bash ${{base_name}}.sh
        '''

rule summary_stairway:
    input:
        expand(f'{workDir}/{vcf_f}/Stairway_results/sampling{{num}}/sampling{{num}}.final.summary',num=range(config["numSampling"]))
    output:
        f'{workDir}/{vcf_f}/Stairway_results.json'
    params:
        generation_interval=config["generation_interval"]
    run:
        import os
        import json
        combined_li={}
        files = input

        for file in files:
            Ne = []  # Store Ne estimates for this file
            with open(file, "r") as f:
                for line in f:
                    if line.startswith("mutation_per_site"):  # Skip header
                        continue
                    values = line.strip().split("\t")
                    try:
                        T = float(values[5])/params.generation_interval  # 'year' column
                        N = float(values[6])  # 'Ne_median' column
                        Ne.append([T, N])
                    except ValueError:
                        print(f"Skipping invalid line in {file}: {line}")          
            combined_li[file]=dict(Ne)
        with open(output[0], 'w') as json_file:
            json.dump(combined_li, json_file, indent=4)






        




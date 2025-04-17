configfile: "./config/Relate_config.yaml"

rule prepare_popLable:
    input:
        f'{workDir}/{vcf_f}/sampling{{num}}.popLabel'
    output:
        f'{workDir}/{vcf_f}/Relate_results/sampling{{num}}.poplabels'
    run:
        import pandas as pd
        import os
        for i in range(len(input)):
            df = pd.read_csv(input[i], sep="\t")  
            # Add new columns
            df['group'] = df['population']
            df['sex'] = "NA"
            df.to_csv(output[i], sep="\t", index=False)

rule prepare_relate_format:
    input:
        f'{workDir}/{vcf_f}/sampling{{num}}.vcf.gz',
    output:
        f'{workDir}/{vcf_f}/Relate_results/sampling{{num}}.sample',
        f'{workDir}/{vcf_f}/Relate_results/sampling{{num}}.haps',
    params:
        relate=config["Relate_path"]
    shell:
        '''
        input_file="{input}"
        input_name=$(basename {input} .vcf.gz)
        work_dir=$(dirname {input})
        output_dir="$work_dir/Relate_results"

        # Create output directory if it doesn't exist
        mkdir -p $output_dir
        cd $output_dir

        # Convert VCF to Relate format
        {params.relate}/bin/RelateFileFormats \
            --mode ConvertFromVcf \
            --haps $input_name.haps \
            --sample $input_name.sample \
            -i $work_dir/$input_name
        '''


rule prepare_relate:
    input:
        haps=f'{workDir}/{vcf_f}/Relate_results/{{sample_file}}.haps',
        sample=f'{workDir}/{vcf_f}/Relate_results/{{sample_file}}.sample',
        poplabels=f'{workDir}/{vcf_f}/Relate_results/{{sample_file}}.poplabels'
    output:
        f'{workDir}/{vcf_f}/Relate_results/{{sample_file}}.annot'
    params:
        relate=config["Relate_path"],
        ancestral = config["ancestral_genome"]
    shell:
        '''
        work_dir=$(dirname {input.haps})
        cd $work_dir
        if [ -f "*.gz" ]; then rm "*.gz"; fi
        # if data is not clear, you need input mask and carry on "remove non-biallelic"
        # {params.relate}/scripts/PrepareInputFiles/PrepareInputFiles.sh
        {params.relate}/bin/RelateFileFormats \
            --mode GenerateSNPAnnotations \
            --haps {input.haps} \
            --sample {input.sample} \
            --ancestor {params.ancestral} \
            --poplabels {input.poplabels} \
            -o {wildcards.sample_file}
        '''

rule run_relate:
    input:
        haps=f'{workDir}/{vcf_f}/Relate_results/{{sample_file}}.haps',
        sample=f'{workDir}/{vcf_f}/Relate_results/{{sample_file}}.sample',
        annot=f'{workDir}/{vcf_f}/Relate_results/{{sample_file}}.annot'
    output:
        f'{workDir}/{vcf_f}/Relate_results/{{sample_file}}.anc',
        f'{workDir}/{vcf_f}/Relate_results/{{sample_file}}.mut'
    params:
        relate=config["Relate_path"],
        genMap=f"{vcf}.genetic_map.txt",
        mu=config["mutation_rate"],
        N=config["N"]
    shell:
        '''
        work_dir=$(dirname {input.haps})
        cd $work_dir
        {params.relate}/scripts/RelateParallel/RelateParallel.sh\
            -m {params.mu} \
            -N  {params.N} \
            --haps {input.haps} \
            --sample {input.sample} \
            --map {params.genMap} \
            --annot {input.annot} \
            --seed 1 \
            --threads 8 \
            -o  {wildcards.sample_file}
            
        '''


rule run_relate_Ne:
    input:
        anc=f'{workDir}/{vcf_f}/Relate_results/{{sample_file}}.anc',
        mut=f'{workDir}/{vcf_f}/Relate_results/{{sample_file}}.mut',
        poplabels=f'{workDir}/{vcf_f}/Relate_results/{{sample_file}}.poplabels'
    output:
        f'{workDir}/{vcf_f}/Relate_results/{{sample_file}}_popsize.coal'
    params:
        relate=config["Relate_path"],
        mu=config["mutation_rate"],
        N=config["N"]
    shell:
        '''
        work_dir=$(dirname {input.anc})
        cd $work_dir
        {params.relate}/scripts/EstimatePopulationSize/EstimatePopulationSize.sh \
            -i {wildcards.sample_file} \
            -m {params.mu} \
            --poplabels {input.poplabels} \
            --seed 1 \
            -o {wildcards.sample_file}_popsize
        '''

rule summary_Relate:
    input:
        expand(f'{workDir}/{vcf_f}/Relate_results/sampling{{sample}}_popsize.coal',sample=range(config["numSampling"]))
    output:
        f'{workDir}/{vcf_f}/Relate_results.json'
    run:
        import pandas as pd
        import numpy as np
        # Define the directory to search
        search_directory = input
        # List to store found files
        combined_li = {}

        for filename in search_directory:
            groups = pd.read_csv(filename, sep="\s+", header=None, nrows=1).values.flatten()
            t = pd.read_csv(filename, sep="\s+", header=None, skiprows=1, nrows=1).values.flatten()
            num_pops = int(np.round(np.sqrt(pd.read_csv(filename, sep="\s+", header=None, skiprows=2).shape[0])))
            results={}
            for p1 in range(1, num_pops + 1):
                for p2 in range(1, p1 + 1):
                    skip_rows = (p1 - 1) * num_pops + p2 + 1
                    c = pd.read_csv(filename, sep="\s+", header=None, skiprows=skip_rows, nrows=1).values.flatten()[2:]
                    group_str = f"{groups[p1 - 1]} - {groups[p2 - 1]}"
                    results[group_str]=dict(zip(t,0.5 / c))

            combined_li[filename]=results

        with open(output[0], 'w') as json_file:
            json.dump(combined_li, json_file, indent=4)

 
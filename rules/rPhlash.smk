# Necessary file:.vcf.gz,poplabel
# option parameters: mutation_rate

rule prepare_relate:
    input:
        f'{workDir}/{vcf_f}/sampling{{num}}.vcf.gz'
    output:
        f'{workDir}/{vcf_f}/Relate_results/sampling{{num}}.sample',
        f'{workDir}/{vcf_f}/Relate_results/sampling{{num}}.haps'
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

rule prepare_popLable:
    input:
        f'{workDir}/{vcf_f}/sampling{{num}}.poplabels'
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

rule prepare_relate:
    input:
        haps=f'{workDir}/{vcf_f}/Relate_results/{{sample_file}}.haps',
        sample=f'{workDir}/{vcf_f}/Relate_results/{{sample_file}}.sample',
        poplabels=f'{workDir}/{vcf_f}/Relate_results/{{sample_file}}.poplabels'
    output:
        f'{workDir}/{vcf_f}/Relate_results/{{sample_file}}_prepared.haps',
        f'{workDir}/{vcf_f}/Relate_results/{{sample_file}}_prepared.sample',
        f'{workDir}/{vcf_f}/Relate_results/{{sample_file}}_prepared.annot'
    params:
        relate=config["Relate_path"],
        ancestal = config["ancestral_genome"]
    shell:
        '''
        work_dir=$(dirname {input})
        cd $work_dir
        {params.relate}/scripts/PrepareInputFiles/PrepareInputFiles.sh \
            --haps {input.haps} \
            --sample {input.sample} \
            --ancestor {params.ancestral} \
            --poplabels {input.poplabels} \
            -o {wildcards.sample_file}_prepared
        # Remove intermediate files
        rm {input.haps}
        rm {input.sample}
        '''

rule run_relate:
    input:
        haps=f'{workDir}/{vcf_f}/Relate_results/{{sample_file}}_prepared.haps',
        sample=f'{workDir}/{vcf_f}/Relate_results/{{sample_file}}_prepared.sample',
        annot=f'{workDir}/{vcf_f}/Relate_results/{{sample_file}}_prepared.annot'
    output:
        f'{workDir}/{vcf_f}/Relate_results/{{sample_file}}.anc',
        f'{workDir}/{vcf_f}/Relate_results/{{sample_file}}.mut'
    params:
        relate=config["Relate_path"],
        genMap={vcf}.genMap.txt,
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
            -o  {wildcards.sample_file} \
            --threads {threads}
        '''


rule run_relate:
    input:
        anc=f'{workDir}/{vcf_f}/Relate_results/{{sample_file}}.anc',
        mut=f'{workDir}/{vcf_f}/Relate_results/{{sample_file}}.mut',
        poplabels=f'{workDir}/{vcf_f}/Relate_results/{{sample_file}}.poplabel'
    output:
        f'{workDir}/{vcf_f}/Relate_results/{{sample_file}}_popsize.coal'
    params:
        relate=config["Relate_path"],
        mu=config["mutation_rate"]
        N=config["N"]
    shell:
        '''
        work_dir=$(dirname {input.anc})
        cd $work_dir
        {params.relate}/scripts/EstimatePopulationSize/EstimatePopulationSize.sh \
            -i {wildcards.sample_file} \
            -m {params.mu} \
            --poplabels {input.poplabel} \
            --seed 1 \
            -o {wildcards.sample_file}_popsize
        '''

#rule summary_tsinfer:
    #input:
        #expand(f'{workDir}/{vcf_f}/tskit_results/sampling{{sample}}.dated.trees',sample=range(config["numSampling"]))
    #output:
        #f'{workDir}/{vcf_f}/Tskit_results.json'
    #params:
        #timepoint = config["timepoint"],
    #conda: "tsinferHan"
    #script:
        "../script/ts_summary_Ne.py"
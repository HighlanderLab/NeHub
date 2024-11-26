configfile: "./config/tskit_config.yaml"

rule run_prepare_ts:
    input:
        f'{workDir}/{vcf_f}/sampling{{num}}.vcf.gz'
    output:
        #expand(f'{workDir}/{vcf_f}/tskit_results/sampling{{number}}.samples', number=range(config["numSampling"]))
        f'{workDir}/{vcf_f}/tskit_results/sampling{{num}}.trees'
    conda: "HLab_tsinfer"
    params:
        ancestal = config["ancestral_genome"]
    shell:
        '''
        parallel -j {threads} '
        python ./script/ts_prepare.py --vcf {{}} --ancestral {params.ancestal}
        ' ::: {input}
        '''
        

rule run_tsinfer:
    input:
        #expand(f'{workDir}/{vcf_f}/tskit_results/sampling{{number}}.samples', number=range(config["numSampling"]))
        f'{workDir}/{vcf_f}/tskit_results/{{sample_file}}.samples'
    output:
        f'{workDir}/{vcf_f}/tskit_results/{{sample_file}}.dated.trees'
    conda: "treeStuff"
    params:
        mu=config["mutation_rate"]
    shell:
        '''
        parallel -j {threads} '
        python ./script/ts_date.py --input {{}} --mutation_rate {params.mu}
        ' ::: {input}
        '''


rule summary_tsinfer:
    input:
        expand(f'{workDir}/{vcf_f}/tskit_results/sampling{{sample}}.dated.trees',sample=range(config["numSampling"]))
    output:
        f'{workDir}/{vcf_f}/Tskit_results.json'
    params:
        timepoint = config["timepoint"],
    conda: "HLab_tsinfer"
    script:
        "../script/ts_summary_Ne.py"
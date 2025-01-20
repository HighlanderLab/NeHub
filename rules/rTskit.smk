configfile: "./config/tskit_config.yaml"

rule run_prepare_ts:
    input:
        vcffile=f'{workDir}/{vcf_f}/sampling{{num}}.vcf.gz'
    output:
        f'{workDir}/{vcf_f}/tskit_results/sampling{{num}}.trees'
    conda: "HLab_tsinfer"
    params:
        ancestal = config["ancestral_genome"]
    shell:
        '''
        python ./script/ts_prepare.py --vcf {input.vcffile} --ancestral {params.ancestal}
        '''
        

rule run_tsinfer:
    input:
        tree=f'{workDir}/{vcf_f}/tskit_results/{{sample_file}}.trees'
    output:
        f'{workDir}/{vcf_f}/tskit_results/{{sample_file}}.dated.trees'
    conda: "tsinferHan"
    params:
        mu=config["mutation_rate"]
    shell:
        '''
        python ./script/ts_date.py --input {input.tree} --mutation_rate {params.mu}
        '''


rule summary_tsinfer:
    input:
        expand(f'{workDir}/{vcf_f}/tskit_results/sampling{{sample}}.dated.trees',sample=range(config["numSampling"]))
    output:
        f'{workDir}/{vcf_f}/Tskit_results.json'
    params:
        timepoint = config["timepoint"],
    conda: "tsinferHan"
    shell:
        '''
        if [ -z {params.timepoint} ]; then
            python ./script/ts_summary_Ne.py --input {input} --output {output[0]}
        else
            python ./script/ts_summary_Ne.py --input {input} --timepoint {params.timepoint} --output {output[0]}
        fi
        '''
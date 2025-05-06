# Necessary file:.vcf.gz,poplabel
# option parameters: mutation_rate

rule prepare_phlash:
    input:
        vcf_phalsh = f'{workDir}/{vcf_f}/sampling{{num}}.vcf.gz',
    output:
        f'{workDir}/{vcf_f}/sampling{{num}}.vcf.gz.tbi'
    conda: "bcftools"
    shell:
        '''
        directory_name=$(dirname {input.vcf_phalsh})
        cd $directory_name
        tabix -f {input.vcf_phalsh}       
        '''

rule run_phlash:
    input:
        vcf_phalsh = f'{workDir}/{vcf_f}/sampling{{num}}.vcf.gz',
        tabix_phlash = f'{workDir}/{vcf_f}/sampling{{num}}.vcf.gz.tbi',
        popLabels = f'{workDir}/{vcf_f}/sampling{{num}}.popLabel'
    output:
        directory(f'{workDir}/{vcf_f}/Phlash_results/sampling{{num}}')
    params:
        mu=config["mutation_rate"]
    conda: 
        'ne_pipe_phlash_nogpu_rocky'
    shell:
        '''
        mkdir -p {output}
        python "./script/phlash_run.py"  --vcf {input.vcf_phalsh} --poplabels {input.popLabels} --out {output}  --mutation_rate {params.mu}
        '''

rule summary_Phlash:
    input:
        expand(f'{workDir}/{vcf_f}/Phlash_results/sampling{{sample}}/',sample=range(config["numSampling"]))
    output:
        f'{workDir}/{vcf_f}/Phlash_results.json'
    conda: 
        'ne_pipe_phlash_nogpu_rocky'
    shell:
        '''
        python "./script/phlash_summary.py"  --input {input} --out {output}
        '''
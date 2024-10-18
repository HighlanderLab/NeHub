rule sampling:
    input:
        f'{vcf}.vcf.gz',
        f'{vcf}.popLabel'
    output:
        expand(f'{workDir}/{vcf_f}/sampling{{number}}.vcf.gz', number=range(config["numSampling"]))
    conda: "bcftools"
    params:
        sampling_num=config["numSampling"],
        sample_size=config["sampleSize"],
        outputDir=lambda wildcards: f"{workDir}/{vcf_f}"
    shell:  
        '''
        if [ ! -d {params.outputDir} ]; then
        mkdir -p "{params.outputDir}"
        fi

        python3 ./script/sampling.py -i {input[1]} -n {params.sampling_num} -s {params.sample_size} -o {params.outputDir}

        # Find all sample*.txt files and run bcftools view on each
        for sample_file in {params.outputDir}/sampling*.txt; do
            base_name=$(basename "$sample_file" .txt)
            output_vcf="{params.outputDir}/${{base_name}}.vcf.gz"
            # Run bcftools view to extract the sample from the VCF
            bcftools view -S "${{sample_file}}" {input[0]} -o $output_vcf
        done
        '''
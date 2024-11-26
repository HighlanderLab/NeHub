rule sampling:
    input:
        f'{vcf}.vcf.gz',
        f'{vcf}.popLabel'
    output:
        expand(f'{workDir}/{vcf_f}/sampling{{number}}.vcf.gz', number=range(config["numSampling"])),
        expand(f'{workDir}/{vcf_f}/sampling{{number}}.popLabel', number=range(config["numSampling"]))
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
        for sample_file in {params.outputDir}/sampling*.popLabel; do
            base_name=$(basename "$sample_file" .popLabel)
            dir_name=$(dirname "$sample_file")
            cd $dir_name
            output_vcf="$base_name".vcf.gz
            awk '{{print $1}}' "$sample_file" > "$base_name".txt
            sed -i '1d' "$base_name".txt
            # Run bcftools view to extract the sample from the VCF
            bcftools view -S "$base_name".txt {input[0]} -Oz -o $output_vcf
            rm "$base_name".txt
        done
        '''
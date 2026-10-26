#!/usr/bin/env nextflow

/*
 * GONE (Genealogy of a population) module
 * 
 * This module implements the GONE method for estimating effective population size (Ne) 
 * using linkage disequilibrium patterns. GONE analyzes genetic data to infer historical 
 * population sizes by examining the decay of linkage disequilibrium over time.
 * 
 * Key features:
 * - Processes VCF files split by population
 * - Handles multiple chromosomes including X and Y
 * - Combine vcf by population and convert to plink format
 * - Generates comprehensive Ne estimates and summary statistics
 * 
 * Input: VCF files split by population
 * Output: 
 *   - Ne estimates (Output_Ne_*)
 *   - d2 estimates (Output_d2_*)
 *   - Detailed output files (OUTPUT_*)
 *   - Summary JSON file (GONE_results.json)
 * 
 * Dependencies:
 *   - bcftools >= 1.9
 *   - plink >= 1.9
 *   - python >= 3.6
 *   - GONE software (https://github.com/esrud/GONE)
 */

workflow run_GONE {
    take:
    pop_chr_files  // Input VCF files grouped by population and chromosome
    
    emit:
    "${params.output_dir}/GONE_results.json"  // Summary results in JSON format
    "${params.output_dir}/GONE_results/Output_d2_*"  // d2 estimates
    "${params.output_dir}/GONE_results/Output_Ne_*"  // Ne estimates
    "${params.output_dir}/GONE_results/OUTPUT_*"  // Detailed output files

    main:
    // Validate GONE parameters to ensure all required settings are present
    if (!params.GONE) {
        error "GONE parameters not defined. Please provide GONE configuration in your params file."
    }
    
    // Define and check required parameters for GONE analysis
    required_params = ['PHASE', 'cMMb', 'DIST', 'NGEN', 'NBIN', 'MAF', 'ZERO', 
                    'maxNCHROM', 'maxNSNP', 'hc', 'REPS', 'threads', 'installation_path']
    missing_params = required_params.findAll { !params.GONE.containsKey(it) }
    if (missing_params) {
        error "Missing required GONE parameters: ${missing_params.join(', ')}"
    }

    // Process and organize VCF files by population and chromosome
    vcfs = pop_chr_files.flatten()
    .map { file -> 
            // Extract population name and chromosome from filename using regex
            def m = file.getName() =~ /^([^_]+_[^_]+)(?:_[^_]+)*_?chr([^\.]+)\.vcf\.gz$/
            if (!m) error "Invalid filename format: ${file.name}. Expected format: population_group.chr*.vcf.gz"
            def population = m[0][1]
            def chr = m[0][2].toUpperCase()  // Convert to uppercase for consistent handling
            tuple(population, chr, file)
        }

    filter_multi_allele(vcfs).set{filtered_vcfs}
    convert2plink(filtered_vcfs).set{filtered_plink}
    run_LD(filtered_plink).set{pop_chr_LD}
    pop_chr_LD.groupTuple(by: 0).map { population, chromosomes, files ->
        [population, chromosomes, files.flatten().collect { file(it) }]
    }.set{chr_LD}
    run_GONE_process(chr_LD).set{GONE_results}
    
    // Verify GONE results to ensure analysis completed successfully
    if (!GONE_results.Ne) {
        error "GONE analysis failed to produce Ne estimates"
    }
    if (!GONE_results.d2) {
        error "GONE analysis failed to produce d2 estimates"
    }
    if (!GONE_results.output) {
        error "GONE analysis failed to produce output files"
    }

    // Generate summary of results
    gone_output_files = GONE_results.Ne.collect()
    summarize_GONE(gone_output_files)
}

process filter_multi_allele {
    tag "${pop}_chr${chr}"
    label "bcftools"
    errorStrategy 'retry'
    maxRetries 1

    input:
    tuple val(pop), val(chr), path(vcf_file)

    output:
    tuple val(pop), val(chr), path("${pop}_chr${chr}_filtered.vcf.gz")

    script:
    """
    ACTIVE_ENV=\$(conda info | grep 'active environment' | awk '{print \$NF}')
    if [ "\$ACTIVE_ENV" = "None" ] || [ "\$ACTIVE_ENV" != "bcftools" ]; then
        source \$(conda info --base)/etc/profile.d/conda.sh
        conda activate bcftools
    fi

    # rename sample IDs: replace _ with X, replace - with x
    bcftools query -l ${vcf_file} > original_samples.txt
    sed 's/_/X/g; s/-/x/g' original_samples.txt > new_samples.txt

    # reheader with new sample names
    bcftools reheader -s new_samples.txt ${vcf_file} -o reheadered.vcf.gz

    # filter multiallelic
    bcftools view -m2 -M2 reheadered.vcf.gz \
        -Oz -o ${pop}_chr${chr}_filtered.vcf.gz
    """
}

process convert2plink {
    tag "${pop}_chr${chr}"
    label "basic_Ne"
    errorStrategy 'retry'
    maxRetries 1
    
    input:
    tuple(val(pop), val(chr),path(vcf_file))

    output:
    tuple(val(pop), val(chr),path("chromosome*.map"),path("chromosome*.ped"))

    script:
    """
    ACTIVE_ENV=\$(conda info | grep 'active environment' | awk '{print \$NF}')
    if [ "\$ACTIVE_ENV" = "None" ] || [ "\$ACTIVE_ENV" != "basic_Ne" ]; then
        source \$(conda info --base)/etc/profile.d/conda.sh
        conda activate basic_Ne
    fi

    plink --vcf ${vcf_file} \
        --make-bed \
        --out "data" \
        --memory 16000 \
        --chr-set 95 \
        --allow-extra-chr

    if [ ! -d "${params.GONE.installation_path}/Linux" ]; then
        echo "Error: GONE installation not found at ${params.GONE.installation_path}/Linux/PROGRAMMES/MANAGE_CHROMOSOMES2"
        exit 1
    fi
    
    echo "run sample variants"
    snp_count=\$(wc -l < data.bim)
    maxNSNP=${params.GONE.maxNSNP}
    if [ \$snp_count -le \$maxNSNP ]; then
    plink --bfile data --export ped --out chromosome${chr}
    else
    prob=\$(awk -v m=\$maxNSNP -v n=\$snp_count 'BEGIN { print m / n }')
    plink --bfile data --thin \$prob -export ped --out chromosome${chr}
    fi

    """
}
/*
 * Process to merge VCF files for each population
 * Uses bcftools to concatenate VCF files while maintaining proper header information
 */
process run_LD {
    tag "${pop}_chr${chr}"
    label "bcftools"
    errorStrategy 'retry'
    maxRetries 2
    
    input:
    tuple(val(pop), val(chr), path(map), path(ped))
    
    output:
    tuple(val(pop), val(chr), path("*"))
    
    script:
    """
    # check GONE installation path
    if [ ! -d "${params.GONE.installation_path}" ]; then
        echo "ERROR: GONE installation_path not found: ${params.GONE.installation_path}"
        echo "Please check your params file or container has GONE installed at this path"
        exit 1
    fi

    # check executable exists
    if [ ! -f "${params.GONE.installation_path}/Linux/PROGRAMMES/LD_SNP_REAL3" ]; then
        echo "ERROR: GONE executable not found at: ${params.GONE.installation_path}/Linux/PROGRAMMES/LD_SNP_REAL3"
        ls ${params.GONE.installation_path} || echo "Directory is empty or inaccessible"
        exit 1
    fi

    SAM=\$(grep -w "" -c ${ped})
    PHASE=${params.GONE.PHASE}
    NGEN=${params.GONE.NGEN}
    NBIN=${params.GONE.NBIN}
    ZERO=${params.GONE.ZERO}
    MAF=${params.GONE.MAF}
    DIST=${params.GONE.DIST}
    cMMb=${params.GONE.cMMb}
    threads=${params.GONE.threads}
    if [ \$threads -eq -99 ]
    then
    threads=\$(getconf _NPROCESSORS_ONLN)
    fi
    options_for_LD="\$SAM \$MAF \$PHASE \$NGEN \$NBIN \$ZERO \$DIST \$cMMb"
    ${params.GONE.installation_path}/Linux/PROGRAMMES/LD_SNP_REAL3 ${chr} \$options_for_LD
    """
}

/*
 * Main GONE analysis process
 * Converts VCF to PLINK format and runs GONE analysis
 * Includes parameter validation and error checking
 */
process run_GONE_process {
    tag "${population}_chr${chromosomes.sort().join('_')}"
    publishDir "${params.output_dir}/GONE", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 1

    input:
    tuple(val(population), val(chromosomes), path(chromosome_files))
    
    output:
    path("Output_d2_*"), emit: d2  // d2 estimates
    path("Output_Ne_*"), emit: Ne  // Ne estimates
    path("OUTPUT_*"), emit: output  // Detailed output
    
    script:
    def sorted_chromosomes = chromosomes.sort()
    def NCHR = sorted_chromosomes.size()
    def FILE = population
    """
    # Check if GONE installation exists
    if [ ! -d "${params.GONE.installation_path}/Linux" ]; then
        echo "Error: GONE installation not found at ${params.GONE.installation_path}/Linux"
        exit 1
    fi
    touch CHROM
    touch OUTPUT
    touch outfileLD
    
    # Process each chromosome
    ${sorted_chromosomes.collect { chr_num ->
        def chr_files = chromosome_files.findAll { it.name.contains(chr_num.toString()) }
        def outfileLD = chr_files.find { it.name.startsWith("outfileLD") }
        def parameters = chr_files.find { it.name.startsWith("parameters") }
        
        """
        # Process chromosome ${chr_num}
        if [ -f "${outfileLD}" ] && [ -f "${parameters}" ]; then
            # Copy files with chromosome-specific naming

            # Append to CHROM file
            cat outfileLD${chr_num} >> CHROM
            
            # Add chromosome header to OUTPUT
            echo "CHROMOSOME ${chr_num}" >> OUTPUT
            
            # Remove lines 2 and 3 from outfileLD
            sed '2,3d' outfileLD${chr_num} > temp
            mv temp outfileLD${chr_num}
            
            # Append parameters to OUTPUT
            cat parameters${chr_num} >> OUTPUT
        fi
        """
    }.join('\n')}

    NGEN=${params.GONE.NGEN}
    NBIN=${params.GONE.NBIN}
    ${params.GONE.installation_path}/Linux/PROGRAMMES/SUMM_REP_CHROM3>>out<<@
    \$NGEN	NGEN
    \$NBIN	NBIN
    ${NCHR}	NCHR
    @

    echo "TOTAL NUMBER OF SNPs" >> OUTPUT_${FILE}
    cat nsnp >> OUTPUT_${FILE}
    echo -e "\\n" >> OUTPUT_${FILE}
    # echo "HARDY-WEINBERG DEVIATION" >> OUTPUT_${FILE}
    # cat outfileHWD >> OUTPUT_${FILE}
    echo -e "\\n" >> OUTPUT_${FILE}
    cat OUTPUT >> OUTPUT_${FILE}
    echo -e "\\n" >> OUTPUT_${FILE}
    echo "INPUT FOR GONE" >> OUTPUT_${FILE}
    echo -e "\\n" >> OUTPUT_${FILE}
    cat outfileLD >> OUTPUT_${FILE}

    rm nsnp
    rm OUTPUT
    rm CHROM

    hc=${params.GONE.hc}
    REPS=${params.GONE.REPS}
    ln -sf "${params.GONE.installation_path}/Linux/PROGRAMMES" ./
    ./PROGRAMMES/GONEparallel.sh -hc \$hc outfileLD \$REPS
    mv outfileLD_Ne_estimates Output_Ne_$FILE
    mv outfileLD_d2_sample Output_d2_$FILE
    # Verify GONE output files
    if [ ! -f "Output_Ne_${FILE}" ]; then
        echo "Error: GONE analysis failed to produce Ne estimates"
        exit 1
    fi
    if [ ! -f "Output_d2_${FILE}" ]; then
        echo "Error: GONE analysis failed to produce d2 estimates"
        exit 1
    fi
    if [ ! -f "OUTPUT_${FILE}" ]; then
        echo "Error: GONE analysis failed to produce output file"
        exit 1
    fi
    """
}

/*
 * Process to summarize GONE results
 * Creates a JSON file with consolidated results from all analyses
 */
process summarize_GONE {
    label "basic_Ne"
    tag "Summarizing GONE results"
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    path(gone_output_files)
    
    output:
    path("GONE_results.json")
    
    script:
    """
    ACTIVE_ENV=\$(conda info | grep 'active environment' | awk '{print \$NF}')
    if [ "\$ACTIVE_ENV" = "None" ] || [ "\$ACTIVE_ENV" != "basic_Ne" ]; then
        source \$(conda info --base)/etc/profile.d/conda.sh
        conda activate basic_Ne
    fi
    
    # Generate summary of results using Python script
    python ${projectDir}/scripts/summarize_results.py \
        -m "GONE" \
        -i ${gone_output_files} \
        -o GONE_results.json

    # Verify summary file was created
    if [ ! -f "GONE_results.json" ]; then
        echo "Error: Failed to create summary file"
        exit 1
    fi
    """
}
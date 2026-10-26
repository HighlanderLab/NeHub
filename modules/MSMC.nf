#!/usr/bin/env nextflow

/*
 * MSMC module for estimating effective population size
 * Uses the Multiple Sequentially Markovian Coalescent (MSMC) method
 * to estimate population size history from genome sequences
 * 
 * This workflow:
 * 1. Converts VCF files to MSMC input format with sampling samples
 * 2. Groups samples by population
 * 3. Runs MSMC analysis for each population
 * 4. Generates a summary of results in JSON format

 * Dependencies:
 * ------------
 * - MSMC2: https://github.com/stschiff/msmc2
 * - Python 3.9+
 * 
 * Input:
 * ------
 * - pop_file: Population definition file (CSV format)
 * - chr_files: Chromosome-specific VCF files
 * 
 * Output:
 * -------
 * - MSMC_results.json: Summary of results in JSON format
 * - MSMC_results/*.final.txt: Detailed MSMC output files
 * 
 * Parameters:
 * ----------
 * Required MSMC parameters in params.MSMC:
 * - sampleSize: Number of samples to process
 * - timeSegmentPattern: Pattern for time segments in analysis
 * - rhoOverMu: Ratio of recombination to mutation rate
 * - maxIterative: Maximum number of iterations
 */

workflow run_MSMC_estimate {
    take:
    pop_file      // Population definition file
    chr_files     // Chromosome-specific VCF files
    
    emit:
    msmc_results

    main:
    // Validate MSMC parameters
    if (!params.MSMC) {
        error "MSMC parameters not defined. Please provide MSMC configuration in your params file."
    }
    // Define required parameters for MSMC analysis
    required_params = ['sampleSize', 'timeSegmentPattern','rhoOverMu','maxIterative']
    // Check if any required parameters are missing from the configuration
    missing_params = required_params.findAll { !params.MSMC.containsKey(it) }
    if (missing_params) {
        error "Missing required MSMC parameters: ${missing_params.join(', ')}"
    }

    // Extract MSMC parameters from configuration
    sampleSize = params.MSMC.sampleSize
    // Convert input files to MSMC format with specified sample size
    chr_file_list = chr_files.collect()
    sampling_convert(pop_file,chr_file_list,sampleSize)
    // Process the output tab files from sampling_convert
    sampling_convert.out.tab_files.collect().flatten()
    .map { file ->
        // Extract population name from filename using regex
        // Example: if filename is "pop1_chr1.tab", extract "pop1"
        def m = file.getName() =~ /^(.*)_chr(.*)_(.*)\.tab$/ 
        if (!m) error "Invalid filename format: ${file.name}"
        // Create tuple of (population_name, file)
        tuple(m[0][1], m[0][3],file)
    }.groupTuple(by: [0, 1])  // Group files by population name
        .set { grouped_tabs }
    // Extract remaining MSMC parameters
    maxIterative = params.MSMC.maxIterative  // Maximum number of iterations for MSMC
    rhoOverMu = params.MSMC.rhoOverMu        // Ratio of recombination to mutation rate
    timeSegmentPattern = params.MSMC.timeSegmentPattern  // Pattern for time segments in analysis

    // Run MSMC analysis on grouped files
    run_MSMC_process(grouped_tabs, maxIterative, rhoOverMu, timeSegmentPattern, sampleSize)
        .collect()  // Collect all results
        .set { msmc_results }

    }

workflow run_MSMC_sum{
    take:
    msmc_results

    main:
    // Generate summary of MSMC results
    summarize_MSMC(msmc_results)
}

workflow run_MSMC{
    take:
    pop_file      // Population definition file
    chr_files     // Chromosome-specific VCF files

    main:
    msmc_results = run_MSMC_estimate(pop_file,chr_files)
    run_MSMC_sum(msmc_results)
    
}

/*
 * Process to convert VCF files to MSMC input format
 * Input: Population file and VCF files
 * Output: Tab-delimited files in MSMC format and CSV labels
 */
process sampling_convert {
    label "basic_Ne"
    tag "Convert VCF to MSMC input"
    publishDir "${params.output_dir}/MSMC", pattern: '*.csv', mode: 'copy'

    input:
    path(pop)         // Population definition file
    path(vcf_file)    // Input VCF file
    val(sampleSize)   // Number of samples to process

    output:
    path("*.tab"), emit: tab_files  // MSMC input files
    path("*.csv"), emit: labels     // Sample labels

    script:
    def seed = params.seed
    """
    ACTIVE_ENV=\$(conda info | grep 'active environment' | awk '{print \$NF}')
    if [ "\$ACTIVE_ENV" = "None" ] || [ "\$ACTIVE_ENV" != "basic_Ne" ]; then
        source \$(conda info --base)/etc/profile.d/conda.sh
        conda activate basic_Ne
    fi

    python ${projectDir}/scripts/vcf2tab_msmc.py \
        -i ${vcf_file} \
        -p ${pop} \
        -s ${sampleSize} \
        -seed ${seed}
    """
}

/*
 * Process to run MSMC analysis
 * Performs demographic inference using the MSMC2 algorithm
 * Input: Grouped tab files by population
 * Output: Final MSMC results for each population
 */
process run_MSMC_process {
    label 'basic_Ne'
    tag "MSMC ${pop}"
    publishDir "${params.output_dir}/MSMC", mode: 'copy'
    input:
    tuple(val(pop), val(seed), path(tab_files))  // Population name and input files
    val(maxIterative)                 // Maximum iterations
    val(rhoOverMu)                    // Recombination/mutation rate ratio
    val(timeSegmentPattern)           // Time segment pattern
    val(sampleSize)                   // Sample size

    output:
    path("${pop}_${seed}.final.txt")          // MSMC results file

    script:
    def input_pattern = (0..<sampleSize*2).join(",")
    """
    ACTIVE_ENV=\$(conda info | grep 'active environment' | awk '{print \$NF}')
    if [ "\$ACTIVE_ENV" = "None" ] || [ "\$ACTIVE_ENV" != "basic_Ne" ]; then
        source \$(conda info --base)/etc/profile.d/conda.sh
        conda activate basic_Ne
    fi

    msmc2_Linux \\
    -I ${input_pattern} \\
    -i ${maxIterative} \\
    -r ${rhoOverMu} \\
    -p ${timeSegmentPattern} \\
    -o ${pop}_${seed} \\
    ${tab_files}
    """
}

/*
 * Process to summarize MSMC results
 * Generates a JSON summary of all MSMC analyses
 * Input: MSMC result files
 * Output: JSON summary file
 */
process summarize_MSMC {
    label "basic_Ne"
    tag "Generate JSON summary"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path msmc_results  // MSMC result files

    output:
    path "MSMC_results.json"  // Summary in JSON format

    script:
    """
    ACTIVE_ENV=\$(conda info | grep 'active environment' | awk '{print \$NF}')
    if [ "\$ACTIVE_ENV" = "None" ] || [ "\$ACTIVE_ENV" != "basic_Ne" ]; then
        source \$(conda info --base)/etc/profile.d/conda.sh
        conda activate basic_Ne
    fi

    python ${projectDir}/scripts/summarize_results.py \\
        -m "MSMC" \\
        -i ${msmc_results} \\
        -o MSMC_results.json \\
        --mutation ${params.mutation_rate}
    """
}
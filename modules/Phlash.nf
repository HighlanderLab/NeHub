#!/usr/bin/env nextflow

/*
 * Phlash module for estimating effective population size
 * 
 * This module implements the PHLASH method to estimate population size history
 * from genetic data. It processes VCF files and population information to
 * generate estimates of effective population size over time.
 * 
 * Output:
 * - JSON files containing population size estimates
 * - Summary results in a consolidated JSON file
 * Dependencies:
 * ------------
 * - PHLASH: Population History Analysis Software
 * - Python 3.9+
 * - Required Python packages:
 *   - numpy: For numerical computations
 *   - pandas: For data manipulation
 *   - scipy: For scientific computing
 *   - matplotlib: For plotting results
 *   - pickle: For saving intermediate results
 * - Required system packages:
 *   - tabix: For indexing and querying VCF files
 *   - bgzip: For compressing VCF files
 * - Required scripts:
 *   - phlash_run.py: Main script for running PHLASH analysis
 *   - summarize_results.py: Generates JSON summary of results
 * 
 * Input requirements:
 * - Population file with sample IDs and population assignments
 * - VCF files containing genetic data
 * - Tabix index files for the VCF files
 * - Mutation rate parameter
 */

workflow run_Phlash {
    take:
    // Input channel containing tuples of (population_file, vcf_file, tabix_file)
    popfile_chr_files

    emit:
    // Main output files
    "${params.output_dir}/Phlash_results.json"  // Consolidated results
    "${params.output_dir}/Phlash_results/*.json"  // Individual population results

    main:
    // Validate required Phlash parameters
    if (!params.Phlash) {
        println "Phlash chr_test not defined, then use the default with 0"
        chr_test=0
    }else{
    required_params = ['chr_test']
    missing_params = required_params.findAll { !params.Phlash.containsKey(it) }
    if (missing_params) {
        println "Missing required Phlash parameters: ${missing_params.join(', ')}, the use the default with 0"
        chr_test=0
    }else{
        chr_test = params.Phlash.chr_test
    }
    }
    // Extract population information from first VCF file
    popfile_chr_files.first().map{
        tuple -> tuple[0]  // Extract population label file from the tuple
    }.set{popfile}

    // Process population labels and group by population
    popfile
        .splitCsv(header: true, sep: '\t')  // Parse tab-separated population file
        .map { row -> 
            // Create unique population_group identifier and map to sample ID
            tuple("${row.population}_${row.group}", row.sample_id) 
        }
        .groupTuple()  // Group samples by their population
        .set { pop_ids }

    // Collect all VCF files for processing
    popfile_chr_files.map{
        tuple -> tuple[1]  // Extract VCF files from the tuple
    }.collect().set{chr_files}

    // Collect all tabix index files
    popfile_chr_files.map{
        tuple -> tuple[2]  // Extract tabix files from the tuple
    }.collect().set{tabix_files}

    // Run Phlash analysis with collected files and mutation rate
    run_Phlash_process(pop_ids, chr_files, tabix_files, params.mutation_rate,chr_test)
        .set { Phalsh_results }

    // Collect all JSON results for summarization
    Phalsh_results.json.collect().set { phlash_json }

    // Generate summary of all results
    summarize_Phlash(phlash_json)
}

// Process to run Phlash analysis for each population
process run_Phlash_process {
    tag "run_Phlash: ${pop_label}"
    label "ne_pipe_phlash"
    publishDir "${params.output_dir}/Phlash", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple(val(pop_label),val(ids))  // Population label and sample IDs
    path(vcf_file)                  // Input VCF files
    path(tabix_file)                // Tabix index files
    val(mutation_rate)              // Mutation rate parameter
    val(chr_test)

    output:
    path("${pop_label}.pkl"), emit: pkl   // Pickle file with raw results
    path("${pop_label}.json"), emit: json // JSON file with processed results

    script:
    """
    ACTIVE_ENV=\$(conda info | grep 'active environment' | awk '{print \$NF}')
    if [ "\$ACTIVE_ENV" = "None" ] || [ "\$ACTIVE_ENV" != "phlash" ]; then
        if ! python -c "import phlash" 2>/dev/null; then
            source \$(conda info --base)/etc/profile.d/conda.sh
            conda activate phlash
        fi
    fi
    export CUDA_VISIBLE_DEVICES=0
    nvidia-smi || echo "No GPU available"
    which python
    python -c "import phlash; print('phlash OK:', phlash.__file__)"
    echo "Processing population: ${pop_label}, group: ${ids}"
    python ${projectDir}/scripts/phlash_run.py \\
        --vcfs ${vcf_file} \\
        --sampleId ${ids.join(' ')} \\
        --out_file ${pop_label} \\
        --chr_test ${chr_test} \\
        --mutation_rate ${mutation_rate}
    """
}

// Process to summarize results from all populations
process summarize_Phlash {
    label "basic_Ne"  // Resource requirements
    publishDir "${params.output_dir}", mode: 'copy'  // Output directory

    input:
    path(Phlash_json)  // Input JSON files from all populations

    output:
    path("Phlash_results.json")  // Consolidated results file

    script:
    """
    ACTIVE_ENV=\$(conda info | grep 'active environment' | awk '{print \$NF}')
    if [ "\$ACTIVE_ENV" = "None" ] || [ "\$ACTIVE_ENV" != "basic_Ne" ]; then
        source \$(conda info --base)/etc/profile.d/conda.sh
        conda activate basic_Ne
    fi
    python ${projectDir}/scripts/summarize_results.py \\
        -m "Phlash" \\
        -i ${Phlash_json} \\
        -o Phlash_results.json  
    """
}


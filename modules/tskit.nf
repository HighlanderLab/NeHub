#!/usr/bin/env nextflow

/*
 * Module: tskit.nf
 * Description: Module for estimating effective population size using tsinfer and tsdate
 * Dependencies:
 *   - Python packages: tsinfer, tsdate, tskit, numpy, pandas
 *   - External tools: Python 3.7+
 * Input:
 *   - VCF files with population labels
 *   - Ancestral genome reference
 * Output:
 *   - Dated tree sequences (.trees.dated)
 *   - Population size estimates (tskit_results.json)

 */

workflow run_tskit {
    take:
    popfile_chr_files  // Input channel containing population-labeled VCF files

    main:
    infer_results = run_tskit_infer(popfile_chr_files)
    run_tskit_ne(infer_results)
}

workflow run_tskit_infer {
    take:
    popfile_chr_files  // Input channel containing population-labeled VCF files

    emit:
    ts_output

    main:
    // Validate tskit parameters
    if (!params.tskit) {
        error "tskit parameters not defined. Please provide tskit configuration in your params file."
    }
    required_params = ['ancestral_genome', 'chunk_size','max_file_size','timepoint']
    missing_params = required_params.findAll { !params.tskit.containsKey(it) }
    if (missing_params) {
        error "Missing required tskit parameters: ${missing_params.join(', ')}"
    }
    if (!file(params.tskit.ancestral_genome).exists()) {
        error "ancestral_genome not found: ${params.tskit.ancestral_genome}"
    }

    // Process input files and extract chromosome information
    popfile_chr_files.map{
        selected_poplabel, vcf, tbx  -> 
        def m= vcf.getName() =~ /^(?:(.+?)[_\.])?chr([^\.]+)\.vcf\.gz$/
        if (!m) error "tskit:: Invalid filename format: ${vcf.name}"
        def prefix = m[0][1] ?: "default"
        tuple(selected_poplabel, vcf, prefix+"_"+m[0][2])
    }.set{files}

    // Set workflow parameters
    ancestral_genome = params.tskit.ancestral_genome
    chunk_size = params.tskit.chunk_size
    max_file_size = params.tskit.max_file_size
    
    // Run tsinfer to create tree sequences
    ts_infer(files,ancestral_genome,chunk_size,max_file_size).set{ts_files}
    
    // Run tsdate to date the tree sequences
    mutation_rate = params.mutation_rate
    ts_date(ts_files,mutation_rate).collect().set{ts_output}
}

workflow run_tskit_ne {
    take:
    dated_trees  // Input channel containing population-labeled VCF files

    main:
    // Generate summary statistics and population size estimates
    summarize_tsinfer(dated_trees)
}

// Process to infer tree sequences from VCF files
process ts_infer {
    tag "tsinfer: ${vcf_file.baseName}"
    label "ts_tools"
    publishDir "${params.output_dir}/tskit", mode: 'copy'

    input:
    tuple(path(poplabel),path(vcf_file),val(output))  // Population label, VCF file, and output prefix
    path(ancestral_genome)  // Reference ancestral genome
    val(chunk_size)  // Size of genomic chunks to process
    val(max_file_size)  // Maximum file size for processing

    output:
    path("${output}.trees.tsz")  // Output tree sequence file

    script:
    """
    ACTIVE_ENV=\$(conda info | grep 'active environment' | awk '{print \$NF}')
    if [ "\$ACTIVE_ENV" = "None" ] || [ "\$ACTIVE_ENV" != "ts_tools" ]; then
        source \$(conda info --base)/etc/profile.d/conda.sh
        conda activate ts_tools
    fi
    python ${projectDir}/scripts/ts_prepare.py --vcf ${vcf_file} --popLabel ${poplabel} \
    --ancestral ${ancestral_genome} --output ${output} \
    --max_file_size ${max_file_size} \
    --chunk_size ${chunk_size} \
    --thread 8
    """
}

// Process to date tree sequences using mutation rate
process ts_date {
    tag "ts_date: ${trees}"
    label "ts_tools"
    publishDir "${params.output_dir}/tskit", mode: 'copy'
    errorStrategy 'ignore'
    
    input:
    path(trees)  // Input tree sequence file
    val(mutation_rate)  // Mutation rate for dating

    output:
    path("*dated.tsz")  // Dated tree sequence file

    script:
    """
    ACTIVE_ENV=\$(conda info | grep 'active environment' | awk '{print \$NF}')
    if [ "\$ACTIVE_ENV" = "None" ] || [ "\$ACTIVE_ENV" != "ts_tools" ]; then
        source \$(conda info --base)/etc/profile.d/conda.sh
        conda activate ts_tools
    fi
    python ${projectDir}/scripts/ts_date.py \
    --mutation_rate ${mutation_rate} \
    --input ${trees} 
    """
}

// Process to summarize results and estimate population sizes
process summarize_tsinfer  {
    tag "ts_summary"
    label "ts_tools"
    publishDir "${params.output_dir}", mode: 'copy'
    errorStrategy 'ignore'

    input:
    path(ts_output)  // Input dated tree sequences

    output:
    path("tskit_results.json")  // JSON file with population size estimates

    script:
    """
    ACTIVE_ENV=\$(conda info | grep 'active environment' | awk '{print \$NF}')
    if [ "\$ACTIVE_ENV" = "None" ] || [ "\$ACTIVE_ENV" != "ts_tools" ]; then
        source \$(conda info --base)/etc/profile.d/conda.sh
        conda activate ts_tools
    fi
    python ${projectDir}/scripts/ts_summary_Ne.py --input ${ts_output} \
    --timepoint ${params.tskit.timepoint} \
    --output tskit_results.json
    """
}

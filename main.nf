#!/usr/bin/env nextflow

/*
 * Ne Estimation Pipeline
 * This pipeline implements various methods for estimating effective population size (Ne)
 * Methods are organized into three categories:
 * Category 1: GONE, SMCpp, stairway (requires sampling and split_pop)
 * Category 2: tskit, PSMC, Relate, Phlash (requires sampling only)
 * Category 3: MSMC (no preprocessing required)
 */

// Input validation
if (!params.input_dir) {
    error "Input directory not specified. Please provide --input_dir parameter"
}

if (!params.output_dir) {
    error "Output directory not specified. Please provide --output_dir parameter"
}

if (!params.chromosomes) {
    error "Chromosomes not specified. Please provide --chromosomes parameter"
}

if (!params.population_file) {
    error "Population file not specified. Please provide --population_file parameter"
}

// Sanitize input parameters
params.input_dir = params.input_dir.toString().replaceAll('/+$', '')
params.output_dir = params.output_dir.toString().replaceAll('/+$', '')

// Define method categories
def cat1_methods = ['GONE', 'SMCpp', 'stairway']
def cat2_methods = ['tskit', 'PSMC', 'Relate', 'Phlash']
def cat3_methods = ['MSMC']

// Parse and validate methods
def methods = (params.method instanceof String)
    ? params.method.tokenize(',').collect { it.trim() }
    : params.method
// set seed
params.seed = params.seed ?: 42
def valid_methods = cat1_methods + cat2_methods + cat3_methods
def invalid_methods = methods.findAll { !valid_methods.contains(it) }
if (invalid_methods) {
    error "Invalid methods specified: ${invalid_methods.join(', ')}. Valid methods are: ${valid_methods.join(', ')}"
}

// Include required modules
if (methods.contains('GONE')) include { run_GONE } from './modules/GONE.nf'
if (methods.contains('SMCpp')) include { run_SMCpp } from './modules/SMCpp.nf'
if (methods.contains('stairway')) include { run_stairway } from './modules/stairway.nf'
if (methods.contains('tskit')) include { run_tskit } from './modules/tskit.nf'
if (methods.contains('tskit')) include { run_tskit_infer } from './modules/tskit.nf'
if (methods.contains('tskit')) include { run_tskit_ne } from './modules/tskit.nf'
if (methods.contains('PSMC')) include { run_PSMC } from './modules/PSMC.nf'
if (methods.contains('Relate')) include { run_Relate } from './modules/Relate.nf'
if (methods.contains('Relate')) include { run_Relate_infer } from './modules/Relate.nf'
if (methods.contains('Relate')) include { run_Relate_ne } from './modules/Relate.nf'
if (methods.contains('Phlash')) include { run_Phlash } from './modules/Phlash.nf'
if (methods.contains('MSMC')) include { run_MSMC } from './modules/MSMC.nf'
if (methods.contains('MSMC')) include { run_MSMC_estimate } from './modules/MSMC.nf'
if (methods.contains('MSMC')) include { run_MSMC_sum } from './modules/MSMC.nf'

// Create output directory if it doesn't exist
if (!file(params.output_dir).exists()) {
    file(params.output_dir).mkdirs()
}
if (!file(params.input_dir).exists()) {
    error "Please Check your input path: ${params.input_dir} "
}

workflow {
    // Partition the requested method list into the three processing categories
    def cat1_selected = methods.findAll { cat1_methods.contains(it) }
    def cat2_selected = methods.findAll { cat2_methods.contains(it) }
    def cat3_selected = methods.findAll { cat3_methods.contains(it) }

    // run_mode controls whether inference and Ne-estimation run together ('full')
    // or in two separate Nextflow invocations ('separate' = inference only,
    // 'sum_run' = Ne-estimation only from previously saved intermediates).
    required_params = ['run_mode']
    missing_params  = required_params.findAll { !params.containsKey(it) }
    if (missing_params) {
        run_mode = "full"
    } else {
        run_mode = params.run_mode
    }

    // ── Full / separate-inference pass ──────────────────────────────────────
    if (run_mode != "sum_run") {

        // Build per-chromosome glob patterns and create the VCF input channel
        chromosome_list = (params.chromosomes instanceof String)
            ? params.chromosomes.tokenize(',')
            : params.chromosomes

        chromosome_patterns = chromosome_list.collect { chr ->
            "${params.input_dir.replaceAll('/+$', '')}/*chr${chr}.vcf.gz"
        }
        input_chromosomes = Channel.fromPath(
            chromosome_patterns,
            checkIfExists: true,
            type: 'file'
        )

        // Validate population file exists before starting any work
        population_file = file("${params.input_dir.replaceAll('/+$', '')}/${params.population_file}")
        if (!population_file.exists()) {
            error "Population file not found: ${population_file}"
        }

        // Downsample the population file if params.sampling < 1 (or > 1 for a
        // fixed count).  The result is a single TSV with the selected samples.
        if (cat1_selected || cat2_selected) {
            def sampling = params.containsKey('sampling') ? params.sampling : 1
            sampling2popfile(population_file, sampling, params.seed).set { selected_poplabel }
        }

        // ── Category 1: GONE, SMCpp, stairway ──────────────────────────────
        // These methods work on population-split VCFs (one file per pop×chr).
        if (cat1_selected) {
            split_pop(input_chromosomes, selected_poplabel).set { pop_chr }

            if (cat1_selected.contains('GONE'))     { run_GONE(pop_chr) }
            if (cat1_selected.contains('SMCpp'))    { run_SMCpp(pop_chr) }
            if (cat1_selected.contains('stairway')) { run_stairway(pop_chr) }
        }

        // ── Category 2: tskit, PSMC, Relate, Phlash ────────────────────────
        // These methods use the full multi-population VCF + a tabix index.
        // The population label is broadcast to every chromosome via .combine().
        if (cat2_selected) {
            selected_poplabel
                .combine(indexVCF(input_chromosomes))
                .set { popfile_chr }

            if (cat2_selected.contains('tskit')) {
                // 'separate' mode: only infer trees now; run Ne step later via sum_run
                if (run_mode == 'separate') { run_tskit_infer(popfile_chr) }
                else                        { run_tskit(popfile_chr) }
            }
            if (cat2_selected.contains('PSMC'))  { run_PSMC(popfile_chr) }
            if (cat2_selected.contains('Relate')) {
                if (run_mode == 'separate') { run_Relate_infer(popfile_chr) }
                else                        { run_Relate(popfile_chr) }
            }
            if (cat2_selected.contains('Phlash')) { run_Phlash(popfile_chr) }
        }

        // ── Category 3: MSMC ────────────────────────────────────────────────
        // MSMC reads the raw VCF directly; no population splitting needed.
        if (cat3_selected) {
            if (run_mode == 'separate') { run_MSMC_estimate(population_file, input_chromosomes) }
            else                        { run_MSMC(population_file, input_chromosomes) }
        }

    // ── Sum-only pass (run_mode == 'sum_run') ────────────────────────────────
    // Re-runs only the Ne-estimation step from previously written intermediates.
    // Useful when the inference jobs succeeded but the summary step failed.
    } else {
        if (methods.contains('tskit')) {
            intermediate_results = Channel
                .fromPath("${params.output_dir}/tskit/*trees.dated.tsz")
                .ifEmpty { error "No dated tree sequences found in ${params.output_dir}/tskit/" }
                .collect()
            run_tskit_ne(intermediate_results)
        }

        if (methods.contains('Relate')) {
            // Reconstruct the chromosome list from the .anc files already on disk
            intermediate_chr = Channel
                .fromPath("${params.output_dir}/Relate/Relate_chr*.anc")
                .map { file ->
                    def m = file.getName() =~ /Relate_chr(.+)\.anc/
                    if (!m) error "Cannot extract chromosome from: ${file}"
                    m[0][1]
                }

            intermediate_file = Channel
                .fromPath("${params.output_dir}/Relate/Relate_chr*")

            def sampling = params.containsKey('sampling') ? params.sampling : 1
            popfile = Channel
                .fromPath("${params.output_dir}/sampling/selected_poplabel${sampling}*")
                .ifEmpty { error "No sampling file found in ${params.output_dir}/sampling/ for sampling=${sampling}" }

            run_Relate_ne(intermediate_chr, intermediate_file, popfile)
        }

        if (methods.contains('MSMC')) {
            intermediate_results = Channel
                .fromPath("${params.output_dir}/MSMC/*.final.txt")
                .ifEmpty { error "No MSMC results found in ${params.output_dir}/MSMC/" }
                .collect()
            run_MSMC_sum(intermediate_results)
        }
    }
}

// Process definitions
process sampling2popfile {
    tag "sampling_${sampling}"
    label "basic_Ne"
    publishDir "${params.output_dir}/sampling", mode: 'copy'

    input:
    path population_file
    val sampling
    val seed
    output:
    path("selected_poplabel${sampling}.txt"), emit: selected_poplabel

    script:
    """
    ACTIVE_ENV=\$(conda info | grep 'active environment' | awk '{print \$NF}')
    if [ "\$ACTIVE_ENV" = "None" ] || [ "\$ACTIVE_ENV" != "basic_Ne" ]; then
        source \$(conda info --base)/etc/profile.d/conda.sh
        conda activate basic_Ne
    fi

    python ${projectDir}/scripts/sampler.py \
        -i ${population_file} \
        -o selected_poplabel${sampling}.txt \
        -s ${sampling} \
        -seed ${seed}
    """
}

process split_pop {
    tag "split_${input_chromosomes.baseName}"
    label "bcftools"
    
    input:
    path input_chromosomes
    path selected_poplabel

    output:
    path("*.vcf.gz"), emit: pop_chr_channel

    script:
    """
    ACTIVE_ENV=\$(conda info | grep 'active environment' | awk '{print \$NF}')
    if [ "\$ACTIVE_ENV" = "None" ] || [ "\$ACTIVE_ENV" != "bcftools" ]; then
        source \$(conda info --base)/etc/profile.d/conda.sh
        conda activate bcftools
    fi

    filename=\$(basename ${input_chromosomes} .vcf.gz)
    
    # Extract unique population-group combinations
    tail -n +2 "${selected_poplabel}" | cut -d '\t' -f2,3 | sort | uniq | while read -r population group; do
        echo "Processing population: \$population, group: \$group"
        
        # Create a group-specific output file
        output_file="\${population}_\${group}.txt"
        
        # Extract samples for this population-group
        awk -v pop="\$population" -v grp="\$group" '\$2 == pop && \$3 == grp {print \$1}' "${selected_poplabel}" > \$output_file
        
        # Process VCF file
        bcftools view -S \$output_file ${input_chromosomes} -Oz -o \${population}_\${group}_\${filename}.vcf.gz
        
        # Clean up temporary file
        rm \$output_file
    done
    """
}

process indexVCF {
    tag "index_${vcf_file.baseName}"
    label "bcftools"
    
    input:
    path vcf_file

    output:
    tuple(path(vcf_file), path("${vcf_file}.*i"))

    script:
    """
    ACTIVE_ENV=\$(conda info | grep 'active environment' | awk '{print \$NF}')
    if [ "\$ACTIVE_ENV" = "None" ] || [ "\$ACTIVE_ENV" != "bcftools" ]; then
        source \$(conda info --base)/etc/profile.d/conda.sh
        conda activate bcftools
    fi

    if [ -e "${params.input_dir}/${vcf_file.getName()}.tbi" ]; then
        cp "${params.input_dir}/${vcf_file.getName()}.tbi" .
    elif [ -e "${params.input_dir}/${vcf_file.getName()}.csi" ]; then
        cp "${params.input_dir}/${vcf_file.getName()}.csi" .
    else
        if ! bcftools index -c ${vcf_file} 2>/dev/null; then
            echo "Tabix failed. Decompressing and re-indexing..."
            gunzip -c ${vcf_file} > ${vcf_file}.unzip
            bgzip -c ${vcf_file}.unzip > ${vcf_file}
            rm ${vcf_file}.unzip
            bcftools index -c ${vcf_file}
        fi
    fi
    """
}

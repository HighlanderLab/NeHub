#!/usr/bin/env nextflow

/*
 * SMCpp module for estimating effective population size
 * 
 * This module implements the SMC++ method to estimate population size history from genomic data.
 * The workflow processes VCF files to estimate historical effective population sizes using SMC++.
 * 
 * Input: VCF files grouped by population and chromosome
 * Output: CSV files with population size estimates and a JSON summary
 * 
 * Dependencies:
 * - SMC++ (via Docker container terhorst/smcpp)
 * - bcftools
 * - Python script for result summarization
 */

workflow run_SMCpp {
    take:
    pop_chr_files
    
    emit:
    "${params.output_dir}/SMCpp_results/*.csv"    // Population size estimates per population
    "${params.output_dir}/SMCpp_results.json"     // Summary of all results in JSON format

    main:
    // Group and sort VCF files by population and chromosome
    vcf_groups = pop_chr_files.flatten()
        .map { file -> 
            // Extract population name and chromosome from filename
            def m = file.getName() =~ /^([^_]+_[^_]+)(?:_[^_]+)*_?chr([^\.]+)\.vcf\.gz$/
            if (!m) error "Invalid filename format: ${file.name}. Expected format: population_group.chr*.vcf.gz"
            def population = m[0][1]
            tuple(population,file)
        }
    
    // Index VCF and extract sample labels (bcftools env)
    index_and_label_smc(vcf_groups).set { indexed_files }

    // Convert indexed VCF to SMC++ format (smcpp env)
    vcf2smc_process(indexed_files).set { smc_files }
    
    // Group SMC files by population for joint analysis
    smc_files
        .map { file ->         
            def m = file.getName()=~ /^(.*)\.chr.*/
            if (!m) error "Invalid filename format: ${file.name}"
            tuple(m[0][1], file)
        }
        .groupTuple()
        .set { grouped_tabs }
    
    // Run SMC++ analysis with specified mutation rate
    run_smcpp_process(grouped_tabs,params.mutation_rate)
        .set { final_json }
    draw_smcpp_process(final_json).collect().set{smcpp_results}
    // Generate summary of all SMC++ results
    summarize_SMCpp(smcpp_results)
}

// Process 1: index VCF and extract sample labels
process index_and_label_smc {
    tag "index and label: ${vcf_file.baseName}"
    stageInMode = 'copy'
    label "bcftools"

    input:
    tuple(val(pop_name), path(vcf_file))

    output:
    tuple(val(pop_name), path(vcf_file), path("*.tbi"), path("*.labels"))

    script:
    """
    ACTIVE_ENV=\$(conda info | grep 'active environment' | awk '{print \$NF}')
    if [ "\$ACTIVE_ENV" = "None" ] || [ "\$ACTIVE_ENV" != "bcftools" ]; then
        source \$(conda info --base)/etc/profile.d/conda.sh
        conda activate bcftools
    fi
    bcftools query -l ${vcf_file} | paste -sd "," > ${pop_name}.labels
    tabix ${vcf_file}
    """
}

// Process 2: convert indexed VCF to SMC++ format
process vcf2smc_process {
    tag "vcf2smc: ${vcf_file.baseName}"
    stageInMode = 'copy'
    label "smcpp"

    input:
    tuple(val(pop_name), path(vcf_file), path(tbi_file), path(labels_file))

    output:
    path("${pop_name}.*.smc.gz")

    script:
    """
    if [ -n "\${SMCPP+x}" ]; then
        echo "Using SMCPP: \$SMCPP"
    elif command -v smc++ &>/dev/null; then
        export SMCPP="smc++"
    else
        export SMCPP="conda run -n smcpp smc++"
    fi
    contig=\$(basename ${vcf_file} | sed -E 's/.*(chr[0-9XY]+).*/\\1/')
    pop_labels=\$(cat ${labels_file})
    \$SMCPP vcf2smc ${vcf_file} ${pop_name}.\${contig}.smc.gz \${contig} pop1:\$pop_labels
    """
}

// Process to run SMC++ analysis
// Estimates population size history using SMC++ method
process run_smcpp_process {
    tag "run smc++ of population: ${pop}"
    errorStrategy 'ignore'
    stageInMode = 'copy'
    label "smcpp"
    
    input:
    tuple(val(pop), path(smc_files))    // Population name and SMC++ format files
    val(mutation_rate)                  // Mutation rate parameter for SMC++

    output:
    tuple(val(pop),path("${pop}/model.final.json"))    // Population size estimates

    script:
    """

    if [ -n "\${SMCPP+x}" ]; then
        echo "Using SMCPP: \$SMCPP"
    elif command -v smc++ &>/dev/null; then
        export SMCPP="smc++"
    else
        export SMCPP="conda run -n smcpp smc++"
    fi
    
    mkdir -p ${pop}
    set +e
    \$SMCPP estimate -o ${pop}/ ${mutation_rate} ${smc_files}
    set -e

    if [[ -s "${pop}/model.final.json" ]]; then
        echo "[INFO] Output detected, continuing pipeline..."
    else
        echo "[WARNING] No valid output detected in ${pop}"
    fi
    """

}
// process run_smcpp_process {
//     tag "run smc++ of population: ${pop}"
//     errorStrategy 'ignore'
//     stageInMode = 'copy'
//     label "smcpp"
    
//     input:
//     tuple(val(pop), path(smc_files))    // Population name and SMC++ format files
//     val(mutation_rate)                  // Mutation rate parameter for SMC++

//     output:
//     tuple(val(pop),path("${pop}/model.final.json"))    // Population size estimates

//     script:
//     """
//     export MPLBACKEND=Agg
//     export MPLCONFIGDIR=/tmp/matplotlib_\$\$
//     mkdir -p \$MPLCONFIGDIR
//     export PYTHONPATH=\${PYTHONPATH:-}

//     if [ -n "\${SMCPP+x}" ]; then
//         echo "Using SMCPP: \$SMCPP"
//     elif command -v smc++ &>/dev/null; then
//         export SMCPP="smc++"
//     else
//         export SMCPP="conda run -n smcpp smc++"
//     fi

//     # monkey patch ascii_plotter via sitecustomize
//     mkdir -p /tmp/sitepkg_\$\$
//     cat > /tmp/sitepkg_\$\$/sitecustomize.py << 'PYEOF'
//     try:
//         import smcpp.optimize.plugins.ascii_plotter as m
//         class _Dummy:
//             def update(self, *args, **kwargs): pass
//             def start(self, *args, **kwargs): pass
//             def finish(self, *args, **kwargs): pass
//         m.AsciiPlotter = _Dummy
//     except Exception:
//         pass
//     PYEOF

//     export PYTHONPATH=/tmp/sitepkg_\$\$:\${PYTHONPATH}

//     mkdir -p ${pop}
//     set +e
//     \$SMCPP estimate -o ${pop}/ ${mutation_rate} ${smc_files}
//     set -e

//     if [[ -s "${pop}/model.final.json" ]]; then
//         echo "[INFO] Output detected, continuing pipeline..."
//     else
//         echo "[WARNING] No valid output detected in ${pop}"
//     fi
//     """

// }

process draw_smcpp_process {
    tag "plot smc++ of population: ${pop}"
    publishDir "${params.output_dir}/SMCpp", mode: 'copy'
    stageInMode = 'copy'
    label "smcpp"
    input:
    tuple(val(pop), path(smc_files))    // Population name and SMC++ format files

    output:
    path("${pop}.csv"), emit: smc_results    // Population size estimates

    script:
    """
    if [ -n "\${SMCPP+x}" ]; then
        echo "Using SMCPP: \$SMCPP"
    elif command -v smc++ &>/dev/null; then
        export SMCPP="smc++"
    else
        export SMCPP="conda run -n smcpp smc++"
    fi

    \$SMCPP plot -c ${pop} model.final.json
    """


}
// Process to summarize SMC++ results
// Combines results from all populations into a single JSON file
process summarize_SMCpp {
    label "basic_Ne"
    tag "Summarizing SMCpp results"
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    path(smcpp_output)    // SMC++ results from all populations
    
    output:
    path("SMCpp_results.json"), emit: summary_json    // Combined summary in JSON format
    
    script:
    """
    ACTIVE_ENV=\$(conda info | grep 'active environment' | awk '{print \$NF}')
    if [ "\$ACTIVE_ENV" = "None" ] || [ "\$ACTIVE_ENV" != "basic_Ne" ]; then
        source \$(conda info --base)/etc/profile.d/conda.sh
        conda activate basic_Ne
    fi

    # Run Python script to summarize results
    python ${projectDir}/scripts/summarize_results.py \
        -m "SMCpp" \
        -i ${smcpp_output} \
        -o SMCpp_results.json

    # Verify summary file was created
    if [ ! -f "SMCpp_results.json" ]; then
        echo "Error: Failed to create summary file"
        exit 1
    fi
    """
}
#!/usr/bin/env nextflow

/*
 * ===========================================================================
 * Module: stairway.nf
 * Description: Stairway Plot 2 module for estimating effective population size
 * 
 * Input:
 *   - pop_chr_files: Population-specific VCF files grouped by chromosome
 * 
 * Output:
 *   - stairway_results.json: Summary of Stairway Plot results
 *   - stairway_results/*.final.summary: Detailed Stairway Plot results per population
 * 
 * Dependencies:
 *   - bcftools: For VCF manipulation and filtering
 *   - vcftools: For VCF thinning
 *   - tabix: For VCF indexing
 *   - bgzip: For VCF compression
 *   - Java: For running Stairway Plot 2
 *   - Python: For running polarize.py and summarize_results.py
 *   - easysfs: For generating site frequency spectra
 * 
 * Required Parameters:
 *   - stairway.easysfs_path: Path to easysfs executable
 *   - stairway.stairway_path: Path to Stairway Plot 2 jar file
 *   - stairway.ancestral_genome: Path to ancestral genome (if fold > 1)
 *   - stairway.fold: Whether to use folded (1) or unfolded (2) SFS
 *   - stairway.ploids: Ploidy level of the samples
 *   - stairway.generation_interval: Years per generation
 *   - stairway.pct_training: Percentage of data for training
 *   - stairway.nrand: Number of randomizations
 *   - stairway.ninput: Number of input files
 *   - stairway.xrange: X-axis range for plotting
 *   - stairway.yrange: Y-axis range for plotting
 *   - stairway.xspacing: X-axis spacing for plotting
 *   - stairway.yspacing: Y-axis spacing for plotting
 *   - stairway.fontsize: Font size for plotting
 *   - mutation_rate: Mutation rate per base pair per generation
 * ===========================================================================
 */

workflow run_stairway {
    take:
    pop_chr_files
    
    emit:
    "${params.output_dir}/stairway_results.json"
    "${params.output_dir}/stairway_results/*.final.summary"

    main:
    // Validate Stairway parameters
    if (!params.stairway) {
        error "stairway parameters not defined. Please provide stairway configuration in your params file."
    }
    required_params = ['easysfs_path','stairway_path','ancestral_genome','fold','ploids','generation_interval','pct_training',
    'nrand','ninput','xrange','yrange','xspacing','yspacing','fontsize']
    missing_params = required_params.findAll { !params.stairway.containsKey(it) }
    if (missing_params) {
        error "Missing required stairway parameters: ${missing_params.join(', ')}"
    }
    if (!file(params.stairway.ancestral_genome).exists()) {
        error "ancestral_genome not found: ${params.stairway.ancestral_genome}"
    }

    // Group by population
    // Extract population name and group VCF files by population
    pop_vcf = pop_chr_files.flatten()
    .map { file -> 
            // Extract population name from filename (e.g., "pop1.chr1.vcf.gz" -> "pop1")
            def m = file.getName() =~ /^([^_]+_[^_]+)(?:_[^_]+)*_?chr([^\.]+)\.vcf\.gz$/
            if (!m) error "Invalid filename format: ${file.name}. Expected format: population_group.chr*.vcf.gz"
            def population = m[0][1]
            tuple(population,file)
        }


    if (params.stairway.fold == 1){
        pop_vcf_files = pop_vcf
    }else{
        if (!params.stairway.ancestral_genome) {
        print("Required parameter 'params.stairway.ancestral_genome' is not defined. use reference genome as ancestral genome directly")
        pop_vcf_files = pop_vcf
        }else{
        polarize(pop_vcf, params.stairway.ancestral_genome).set{pop_vcf_files}
        
        }
    }
    def total_length = params.stairway.total_length ?: "NA"
    filter_missing(pop_vcf_files,total_length).set{filtered_vcf}
    generate_sfs(filtered_vcf).set{sfs}
    sfs.groupTuple(by: 0).map { pop, poplabel, lengths, sfs_lists ->
        def total_len = lengths.sum()
        tuple(pop, poplabel.first(), total_len, sfs_lists.flatten())
    }.set{sfs_groups}
    // Group VCF files by population
    group_by_population_sfs(sfs_groups,total_length).set{combined_sfs}

    write_stairway_blueprints(combined_sfs).set{blueprint}
    
    run_stairway_process(blueprint)
    .collect()
    .ifEmpty { error "Stairway analysis failed to produce results" }
    .set{stairway_results}
    summarize_stairway_results(stairway_results)
}

/*
 * Process: polarize
 * Purpose: Polarizes VCF files using ancestral genome information
 * Input:
 *   - pop: Population identifier
 *   - vcf_file: Input VCF file
 *   - ancestral_genome: Reference ancestral genome
 * Output:
 *   - polarized VCF file
 */
process polarize {
    tag "polarize: ${vcf_file.baseName}"
    label "ts_tools"
    input:
    tuple(val(pop),path(vcf_file))
    path(ancestral_genome)

    output:
    tuple(val(pop),path("polarized*"))

    script:
    """
    ACTIVE_ENV=\$(conda info | grep 'active environment' | awk '{print \$NF}')
    if [ "\$ACTIVE_ENV" = "None" ] || [ "\$ACTIVE_ENV" != "ts_tools" ]; then
        source \$(conda info --base)/etc/profile.d/conda.sh
        conda activate ts_tools
    fi
    python ${projectDir}/scripts/polarize.py -v ${vcf_file} -a ${ancestral_genome} -o polarized_${vcf_file}
    """
}

/*
 * Process: filter_missing
 * Purpose: Filters out sites with missing data from VCF files
 * Input:
 *   - pop: Population identifier
 *   - vcf_file: Input VCF file
 * Output:
 *   - Filtered VCF file with missing sites removed
 */
process filter_missing{
    tag "filter missing: ${vcf_file.baseName}"
    label "bcftools"
    input:
    tuple(val(pop),path(vcf_file))
    val(total_length)

    output:
    tuple(val(pop),env(length),path("filtered_${vcf_file}"),path("filtered_${pop}"))

    script:
    """
    ACTIVE_ENV=\$(conda info | grep 'active environment' | awk '{print \$NF}')
    if [ "\$ACTIVE_ENV" = "None" ] || [ "\$ACTIVE_ENV" != "bcftools" ]; then
        source \$(conda info --base)/etc/profile.d/conda.sh
        conda activate bcftools
    fi
    bcftools view -e 'FILTER="MISSING"' ${vcf_file} -Oz -o filtered_${vcf_file}
    bcftools query -l ${vcf_file} | awk '{print \$0, "pop0"}' > "filtered_${pop}"
    if [ "${total_length}" = "NA" ]; then
        length=\$(bcftools view -h ${vcf_file} | \
        grep "^##contig" | \
        awk -F'[,=]' '{print \$5}' | \
        sed 's/>//' | \
        paste -sd+ | \
        bc)
    else
        length=0
    fi
    """
}

/*
 * Process: generate_sfs
 * Purpose: Generates site frequency spectrum from VCF files
 * Input:
 *   - vcf_file: Input VCF file
 *   - poplabel: Population labels file
 * Output:
 *   - Site frequency spectrum files in fastsimcoal2 format
 */
process generate_sfs {
    tag "generate vcf to sfs: ${vcf_file.baseName}"
    label "basic_Ne"
    input:
    tuple(val(pop),val(length),path(vcf_file),path(poplabel))

    output:
    tuple(val(pop),path(poplabel),val(length),path("${vcf_file.getSimpleName()}_SFS.obs"))

    script:
    """
    ACTIVE_ENV=\$(conda info | grep 'active environment' | awk '{print \$NF}')
    if [ "\$ACTIVE_ENV" = "None" ] || [ "\$ACTIVE_ENV" != "basic_Ne" ]; then
        source \$(conda info --base)/etc/profile.d/conda.sh
        conda activate basic_Ne
    fi

    if [ ! -f "${params.stairway.easysfs_path}" ]; then
        echo "ERROR: easysfs_path not found: ${params.stairway.easysfs_path}"
        exit 1
    fi
    
    min_count=\$(wc -l ${poplabel} | awk '{{print \$1}}')
    if [ ${params.stairway.fold} -eq 1 ]; then
        ${params.stairway.easysfs_path} -i ${vcf_file} -p ${poplabel} -a --proj \$min_count -f -o ${poplabel.getName()}_SFS
    else
        project=\$((${params.stairway.ploids} * \${min_count}))
        ${params.stairway.easysfs_path} -i ${vcf_file} -p ${poplabel} -a --proj \$project -o ${poplabel.getName()}_SFS -f --unfolded
    fi
    mv ${poplabel.getName()}_SFS/fastsimcoal2/pop0_DAFpop0.obs ${vcf_file.getSimpleName()}_SFS.obs
    """
}

/*
 * Process: group_by_population_vcf
 * Purpose: Combines and processes VCF files for each population
 * Steps:
 *   1. Concatenates VCF files from different chromosomes
 *   2. Indexes the concatenated VCF
 *   3. Thins the VCF to specified base pair interval
 *   4. Creates population labels file
 * Input:
 *   - pop: Population identifier
 *   - vcfs: List of VCF files to combine
 * Output:
 *   - Combined and thinned VCF file
 *   - Population labels file
 */
process group_by_population_sfs {
    tag "Grouping: ${pop}"
    label "basic_Ne"
    
    input:
    tuple(val(pop),path(poplabel),val(sum_length),path(chr_SFS))
    val(total_length)

    output:
    tuple(env(final_total_length), path(poplabel), path("merged_all_chr.obs"))
    
    script:
    def SFS_list = chr_SFS.collect { "'${it}'" }.join(" ")
    """
    ACTIVE_ENV=\$(conda info | grep 'active environment' | awk '{print \$NF}')
    if [ "\$ACTIVE_ENV" = "None" ] || [ "\$ACTIVE_ENV" != "basic_Ne" ]; then
        source \$(conda info --base)/etc/profile.d/conda.sh
        conda activate basic_Ne
    fi

    if [ "${total_length}" = "NA" ]; then
        final_total_length=${sum_length}
    else
        final_total_length=${total_length}
    fi

    python ${projectDir}/scripts/merge_sfs.py ${SFS_list} -o merged_all_chr.obs
    
    """
}


/*
 * Process: write_stairway_blueprints
 * Purpose: Creates blueprint files for Stairway Plot 2 analysis
 * Steps:
 *   1. Calculates total sequence length
 *   2. Determines number of individuals
 *   3. Processes SFS data
 *   4. Generates blueprint file with all required parameters
 * Input:
 *   - vcf_file: Input VCF file
 *   - poplabel: Population labels file
 *   - SFS_file: Site frequency spectrum file
 * Output:
 *   - Blueprint file for Stairway Plot 2
 */
process write_stairway_blueprints {
    tag "write stairway blueprints of population: ${poplabel}"
    label "bcftools"
    input:
    tuple(val(final_total_length),path(poplabel),path(SFS_file))

    output:
    path("${poplabel}.blueprint")

    script:
    """
    ACTIVE_ENV=\$(conda info | grep 'active environment' | awk '{print \$NF}')
    if [ "\$ACTIVE_ENV" = "None" ] || [ "\$ACTIVE_ENV" != "bcftools" ]; then
        source \$(conda info --base)/etc/profile.d/conda.sh
        conda activate bcftools
    fi
    nind=\$(wc -l < ${poplabel})
    # Calculate sequence count
    nseq=\$(( ${params.stairway.ploids} * \${nind} ))
    if [ ${params.stairway.fold} -eq 1 ]; then
        line3=\$(sed -n '3p' "${SFS_file}" | tr ' ' '\t'|cut -f2-)
    else
        line3=\$(sed -n '3p' "${SFS_file}" | awk '{for(i=2;i<NF;i++) printf "%s%s", \$i, (i<NF-1 ? OFS : ""); print ""}')
    fi
    
    fold=\$(if [ ${params.stairway.fold} -eq 1 ];then echo "true"; else echo "false"; fi)

    # Create parameter file
    cat > ${poplabel}.blueprint << EOL
    popid: ${poplabel}
    nseq: \$nseq
    L: ${final_total_length}
    whether_folded: \$fold
    SFS: \$line3
    pct_training: ${params.stairway.pct_training}
    nrand: ${params.stairway.nrand}
    project_dir:${poplabel}
    stairway_plot_dir: ${params.stairway.stairway_path}
    ninput: ${params.stairway.ninput}
    mu:  ${params.mutation_rate}
    year_per_generation: ${params.stairway.generation_interval}
    plot_title: ${poplabel}
    xrange: ${params.stairway.xrange}
    yrange: ${params.stairway.yrange}
    xspacing: ${params.stairway.xspacing}
    yspacing: ${params.stairway.yspacing}
    fontsize: ${params.stairway.fontsize}
    EOL
    """
}

/*
 * Process: run_stairway_process
 * Purpose: Executes Stairway Plot 2 analysis
 * Steps:
 *   1. Builds Stairway Plot 2 project
 *   2. Runs the analysis
 *   3. Moves results to appropriate location
 * Input:
 *   - blueprint_file: Stairway Plot 2 blueprint file
 * Output:
 *   - Final summary files with population size estimates
 */
process run_stairway_process {
    tag "run stairway: ${blueprint_file.baseName}"
    publishDir "${params.output_dir}/stairway", mode: 'copy'

    input:
    path(blueprint_file)

    output:
    path("*.final.summary")

    script:
    """
    if [ ! -d "${params.stairway.stairway_path}" ]; then
        echo "ERROR: stairway_path not found: ${params.stairway.stairway_path}"
        exit 1
    fi

    directory_name=\$(dirname ${blueprint_file})
    base_name=\$(basename  ${blueprint_file})
    cd \$directory_name
    java -cp ${params.stairway.stairway_path} Stairbuilder  ${blueprint_file}       
    chmod +x \${base_name}.sh
    bash \${base_name}.sh
    folder_path=\$(basename \${base_name} .blueprint)
    mv \${folder_path}/*.final.summary \${folder_path}.final.summary
    """
}

/*
 * Process: summarize_stairway_results
 * Purpose: Aggregates and summarizes Stairway Plot 2 results
 * Input:
 *   - stairway_results: Individual population results
 * Output:
 *   - JSON file containing summarized results
 */
process summarize_stairway_results {
    label "basic_Ne"
    tag "Summarizing stairway results"
    publishDir "${params.output_dir}", mode: 'copy'
    input:
    path(stairway_results)
    
    output:
    path("stairway_results.json")
    
    script:
    """
    ACTIVE_ENV=\$(conda info | grep 'active environment' | awk '{print \$NF}')
    if [ "\$ACTIVE_ENV" = "None" ] || [ "\$ACTIVE_ENV" != "basic_Ne" ]; then
        source \$(conda info --base)/etc/profile.d/conda.sh
        conda activate basic_Ne
    fi

    python ${projectDir}/scripts/summarize_results.py \
        -m "stairway" \
        -i ${stairway_results} \
        --generation_interval ${params.stairway.generation_interval} \
        -o stairway_results.json

    # Verify summary file
    if [ ! -f "stairway_results.json" ]; then
        echo "Error: Failed to create summary file"
        exit 1
    fi
    """
}
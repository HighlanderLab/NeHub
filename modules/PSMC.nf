#!/usr/bin/env nextflow

/*
 * PSMC module for estimating effective population size
 * Uses the Pairwise Sequentially Markovian Coalescent (PSMC) method
 * to estimate population size history from diploid genome sequences
 * 
 * This workflow:
 * 1. Takes population and chromosome files as input
 * 2. Converts VCF files to FASTQ format for PSMC analysis
 * 3. Runs PSMC analysis on each sample
 * 4. Summarizes results into a JSON file
 * 
 * Required parameters in params.PSMC:
 * - ref_genome: Reference genome file path
 * - timeSegmentPattern: PSMC time segment pattern (e.g., "4+25*2+4+6")
 * - bin_size: Size of time bins for analysis
 * 
 * Dependencies:
 * - bcftools: For VCF file manipulation
 * - vcf-consensus: For consensus sequence generation
 * - bgzip: For file compression
 * - psmc: For PSMC analysis
 * - fq2psmcfa: For converting FASTQ to PSMC format
 * - Python: For result summarization
 */

workflow run_PSMC {
    take:
    popfile_chr_files  // Input channel containing population file and chromosome files
    
    emit:
    "${params.output_dir}/psmc_results/*.psmc"  // PSMC result files
    "${params.output_dir}/psmc_results.json"    // Summary JSON file

    main:
    // Validate PSMC parameters
    if (!params.PSMC) {
        error "PSMC parameters not defined. Please provide PSMC configuration in your params file."
    }
    
    // Check for required parameters
    required_params = ['ref_genome', 'timeSegmentPattern','bin_size']
    missing_params = required_params.findAll { !params.PSMC.containsKey(it) }
    if (missing_params) {
        error "Missing required PSMC parameters: ${missing_params.join(', ')}"
    }
    if (!file(params.PSMC.ref_genome).exists()) {
        error "ref_genome not found: ${params.PSMC.ref_genome}"
    }

    // Extract population information from input file
    popfile_chr_files.first().map{
        tuple -> tuple[0]
    }.set{popfile}
    popfile
        .splitCsv(header: true ,sep: '\t')
        .map { row -> tuple("${row.population}_${row.group}", row.sample_id) }
        .set { pop_id }

    // Extract chromosome files
    popfile_chr_files.map{
        tuple -> tuple[1,2]
    }.set{chr_files}

    // Combine population IDs with chromosome files and extract chromosome numbers
    pop_id
        .combine(chr_files)
        .map { pop, id, vcf, tbx -> 
        def m= vcf.getName() =~ /^(?:(.+?)[_\.])?chr([^_\.]+)\.vcf\.gz$/ 
        if (!m) error "PSMC:: Invalid filename format: ${vcf.name}"
        chr=m[0][2]
        tuple(vcf, tbx, id, chr)
        }.set {id_chr}

    // Process files through PSMC pipeline
    ref_genome = params.PSMC.ref_genome
    convert2fq(id_chr,ref_genome).set{fq_files} // Convert VCF to FASTQ format
    
    // Group FASTQ files by sample ID
    fq_files.
        groupTuple()
        .set{id_fq_files}
    mergeFastq(id_fq_files).set{id_fq}
    
    // Run PSMC analysis
    timeSegmentPattern = params.PSMC.timeSegmentPattern
    PSMC(id_fq,timeSegmentPattern ).set { psmc_results }
    
    psmc_results.collect().set{psmc_col}
    // Generate summary of PSMC results
    summarize_PSMC(psmc_col,popfile,params.mutation_rate,params.PSMC.bin_size)
}

/*
 * Process to convert VCF files to FASTQ format for PSMC analysis
 * Input: VCF file, tabix index, sample ID, chromosome number, and reference genome
 * Output: Compressed FASTQ file for the sample and chromosome
 */
process convert2fq {
    tag "vcf to fa (id: ${sample} chr: ${chr})"
    label 'bcftools'
    input:
    tuple(path(vcf_file),path(tab_file),val(sample),val(chr))
    path(ref_genome)
    
    output:
    tuple(val(sample),path("${sample}.chr${chr}.fq.gz"))

    script:
    """
    set -eo pipefail
    ACTIVE_ENV=\$(conda info | grep 'active environment' | awk '{print \$NF}')
    if [ "\$ACTIVE_ENV" = "None" ] || [ "\$ACTIVE_ENV" != "bcftools" ]; then
        source \$(conda info --base)/etc/profile.d/conda.sh
        conda activate bcftools
    fi
    
    # Prepare reference genome
    # Check if file is gzipped
    if [[ ${ref_genome} == *.gz ]]; then
        zcat ${ref_genome} > temp_ref.fa
    else
        cat ${ref_genome} > temp_ref.fa
    fi

    awk -v target="${chr}" '
    /^>/ {
        split(\$1, a, ":");
        sub(/^>/, "", a[1]);
        sub(/^chr/, "", a[1]);
        if (a[1] == target) {
            print ">chr" a[1];
            flag=1
        } else {
            flag=0
        }
        next
    }
    flag { print }
    ' temp_ref.fa > simplified_ref.fa
    rm temp_ref.fa

    cat simplified_ref.fa | bcftools consensus -I -s ${sample} ${vcf_file} | \
            sed 's/>/@/g' | bgzip > ${sample}.chr${chr}.fq.gz
    rm simplified_ref.fa
    """
}

/*
 * Process to merge FASTQ files from different chromosomes for each sample
 * Input: Sample ID and list of chromosome-specific FASTQ files
 * Output: Single merged FASTQ file per sample
 */
process mergeFastq {
    tag "merge Fastq by id"
    label 'bcftools'
    input:
    tuple(val(sample),path(fq_files))  // Collects all per-chr fq.gz files

    output:
    tuple(val(sample),path("${sample}.fq.gz"))  // Final merged file for PSMC

    script:
    """
    set -eo pipefail
    ACTIVE_ENV=\$(conda info | grep 'active environment' | awk '{print \$NF}')
    if [ "\$ACTIVE_ENV" = "None" ] || [ "\$ACTIVE_ENV" != "bcftools" ]; then
        source \$(conda info --base)/etc/profile.d/conda.sh
        conda activate bcftools
    fi

    zcat ${fq_files} | bgzip > ${sample}.fq.gz
    """
}

/*
 * Process to run PSMC analysis
 * Input: Sample ID and merged FASTQ file
 * Output: PSMC result file
 * Parameters:
 * - N25: Number of iterations
 * - t20: Maximum 2N0 coalescent time
 * - r5: Initial theta/rho ratio
 * - p: Time segment pattern
 */
process PSMC {
    tag "run PSMC: ${sample}"
    label 'basic_Ne'
    publishDir "${params.output_dir}/PSMC", mode: 'copy'
    input:
    tuple(val(sample),path(fq_file))
    val(timeSegmentPattern)

    output:
    path("${sample}.psmc")

    script:
    """
    set -eo pipefail
    ACTIVE_ENV=\$(conda info | grep 'active environment' | awk '{print \$NF}')
    if [ "\$ACTIVE_ENV" = "None" ] || [ "\$ACTIVE_ENV" != "basic_Ne" ]; then
        source \$(conda info --base)/etc/profile.d/conda.sh
        conda activate basic_Ne
    fi
    fq2psmcfa -q0 ${fq_file} > ${sample}.psmcfa
    psmc -N25 -t20 -r5 -p ${timeSegmentPattern} ${sample}.psmcfa -o ${sample}.psmc
    """
}

/*
 * Process to summarize PSMC results
 * Input: PSMC result files, group information, mutation rate, and bin size
 * Output: JSON file containing summarized results
 */
process summarize_PSMC {
    tag "Summarizing PSMC results"
    label "basic_Ne"
    publishDir "${params.output_dir}", mode: 'copy'
    input:
    path(psmc_results)
    val(group_file)
    val(mutation_rate)
    val(bin_size)

    output:
    path("PSMC_results.json")

    script:
    """
    ACTIVE_ENV=\$(conda info | grep 'active environment' | awk '{print \$NF}')
    if [ "\$ACTIVE_ENV" = "None" ] || [ "\$ACTIVE_ENV" != "basic_Ne" ]; then
        source \$(conda info --base)/etc/profile.d/conda.sh
        conda activate basic_Ne
    fi
    python ${projectDir}/scripts/summarize_results.py \
        -m "PSMC" \
        -i ${psmc_results} \
        -g ${group_file} \
        --mutation ${mutation_rate} \
        --bin_size ${bin_size} \
        -o PSMC_results.json

    # Verify summary file
    if [ ! -f "PSMC_results.json" ]; then
        echo "Error: Failed to create summary file"
        exit 1
    fi
    """
}
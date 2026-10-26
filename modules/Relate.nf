#!/usr/bin/env nextflow

/*
 * Relate module for estimating effective population size
 * 
 * This module implements the Relate method to estimate population size history
 * from genetic data. It processes VCF files, converts them to Relate format,
 * generates trees, and estimates effective population size (Ne) over time.
 * 
 * Input:
 * - popfile_chr_files: Channel containing population label files and VCF files
 * - ancestral_genome: Ancestral genome reference
 * - genMap_dir: genetic map
 * 
 * Output:
 * - Ancestral tree files (*.anc)
 * - Mutation files (*.mut)
 * - Coalescence files (*.coal)
 * - JSON summary of results
 * 
 * Dependencies:
 * - Relate software
 * - Python 3.x
 * - VCF tools
 */

workflow run_Relate {
    take:
    popfile_chr_files  // Input channel containing population files and VCF files

    main:
    infer_results = run_Relate_infer(popfile_chr_files)
    infer_results.relate_file.flatten().set{all_files}
    run_Relate_ne(infer_results.chr,all_files, infer_results.popfile)

}

workflow run_Relate_ne {
    take:
    chr
    relate_tree_files // Input channel containing population files and VCF files
    popfile

    main:
    if (!params.Relate) {
        error "Relate parameters not defined. Please provide Relate configuration in your params file."
    }
    required_params = ['Relate_path','initial_N']
    missing_params = required_params.findAll { !params.Relate.containsKey(it) }
    if (missing_params) {
        error "Missing required Relate parameters: ${missing_params.join(', ')}"
    }
    output_prefix = "Relate_Ne"
    Relate_path = params.Relate.Relate_path
    relate_tree_files.map { file ->
        def m = file.getName() =~ /Relate_chr([^\.]+)\./
        if (!m) error "Invalid filename format: ${file.name}"
        def chr = m[0][1]
        tuple(chr, file)
    }
    .groupTuple(by:0)    
    .map { chr, files ->
        // Check if any files end with .anc or .mut
        def anc_files = files.findAll { it.getName().endsWith('.anc') }
        def mut_files = files.findAll { it.getName().endsWith('.mut') }
        
        // Add null checks to prevent errors
        if (anc_files.isEmpty()) {
            error "No .anc files found for chromosome ${chr}"
        }
        if (mut_files.isEmpty()) {
            error "No .mut files found for chromosome ${chr}"
        }
        
        tuple(chr, anc_files[0], mut_files[0])
    }
    .set { chr_anc_mut_channel }
    // Remove trees with few mutations for each chromosome
    def threshold = params.Relate.threshold ?: 0.5
    RemoveTrees(Relate_path,output_prefix,threshold,chr_anc_mut_channel)
    chr.collect().set { relate_chr }
    initial_N = params.Relate.initial_N

    def num_iter = params.Relate.num_iter ?: 5
    def num_chunk = params.Relate.num_chunk ?: 10
    // Initialize variables outside the loop
    initial_anc_mut = RemoveTrees.out.anc_mut.collect{chr,anc,mut -> [anc,mut]}
    initial_dist = RemoveTrees.out.dist.collect{chr,dist -> dist}
    CoalRateForTree(Relate_path,output_prefix,relate_chr, initial_anc_mut, initial_dist)

    initial_coal = CoalRateForTree.out.coal
    SampleBranchCallCoal(
        Relate_path,
        output_prefix, 
        num_iter,
        num_chunk,
        params.mutation_rate,
        initial_anc_mut,
        initial_dist,
        initial_coal,
        CoalRateForTree.out.chr_list)

    SampleBranchCallCoal.out.anc.flatMap{it}.set{anc_files}
    anc_files.map{file ->
        def m = file.getName() =~ /_final_chr([^\.]+)\./
        return tuple(m[0][1],file)
    }.set{final_anc}
    SampleBranchCallCoal.out.mut.flatMap{it}.set{mut_files}
    mut_files.map{file ->
        def m = file.getName() =~ /_final_chr([^\.]+)\./
        return tuple(m[0][1],file)
    }.set{final_mut}
    final_anc.combine(final_mut, by:0).set{final_anc_mut}
    SampleBranchCallCoal.out.coal.flatMap{it}.set{coal_files}
    RemoveTrees.out.dist.combine(coal_files).set { final_coal_dist }
    final_anc_mut.combine(final_coal_dist, by: 0).set { pairwise_input }
    output_prefix_final = "Relate_Ne_final"
    CoalRatePairwise(Relate_path,output_prefix_final,CoalRateForTree.out.chr_list,pairwise_input)
    CoalRatePairwise.out.collect().set{bin_file}
    SummarizeCoalRate(Relate_path,output_prefix_final,CoalRateForTree.out.chr_list,bin_file,popfile)
    
    DivideAncMutFinal(Relate_path, output_prefix_final, num_chunk, final_anc_mut)
    DivideAncMutFinal.out.anc_mut.flatMap { chr, files ->
        files.collect { file ->
            def m = file.getName() =~ /_tmp_chr([^\.]+)\./
            if (!m) error "Invalid filename format: ${file.name}"
            def chunk = m[0][1]
            tuple(chr, chunk, file)
        }
    }.groupTuple(by: [0,1]).combine(final_coal_dist, by:0).set{chr_chunk_files_final}
    ReEstimateBranchLengths(Relate_path,output_prefix_final,params.mutation_rate,chr_chunk_files_final )
    ReEstimateBranchLengths.out.groupTuple(by: 0).map{chr,files -> tuple(chr,files.flatten())}.combine(DivideAncMutFinal.out.params_file, by:0).set{chr_files_final}
    CombineAncMutFinal(Relate_path,output_prefix_final,chr_files_final)
    CombineAncMutFinal.out.view()

    // // // Generate summary of results
    summarize_Relate(SummarizeCoalRate.out)
}
    
workflow run_Relate_infer {
    take:
    popfile_chr_files  // Input channel containing population files and VCF files

    emit:
    chr
    relate_file
    popfile

    main:
     // Validate required Relate parameters
    if (!params.Relate) {
        error "Relate parameters not defined. Please provide Relate configuration in your params file."
    }
    required_params = ['ancestral_genome', 'Relate_path','genMap_dir','initial_N']
    missing_params = required_params.findAll { !params.Relate.containsKey(it) }
    if (missing_params) {
        error "Missing required Relate parameters: ${missing_params.join(', ')}"
    }
    if (!file(params.Relate.ancestral_genome).exists()) {
        error "ancestral_genome not found: ${params.Relate.ancestral_genome}"
    }
    if (!file(params.Relate.genMap_dir).exists()) {
        error "genMap_dir not found: ${params.Relate.genMap_dir}"
    }

    // Extract population information from first VCF file
    popfile_chr_files.first().map{
        tuple -> tuple[0]  // Get population label file
    }.set{popfile}

    split_ancestral_chr(params.Relate.ancestral_genome).flatten().map{
        file ->
        def m = file.getName() =~ /(?i)chr([^\.]+)\.fa$/  // Removes .fa extension
        if (!m) error "Invalid filename format: ${file.name}"
        tuple(m[0][1], file)
    }.set { chr_fasta_ch }
    // Process VCF files and extract chromosome information
    popfile_chr_files.map{
        selected_poplabel, vcf, tbx  -> 
        def m = vcf.getName() =~ /^(?:(.+?)[_\.])?chr([^\.]+)\.vcf\.gz$/
        if (!m) error "Relate:: Invalid filename format: ${vcf.name}"
        def prefix = m[0][1] ?: "default"
        tuple(m[0][2],selected_poplabel, vcf, prefix+"_"+m[0][2])
    }.combine(chr_fasta_ch,by:0).set{files}
    // Prepare files for Relate analysis
    Relate_path = params.Relate.Relate_path
    prepare_relate_format(files, Relate_path).set{relate_files}

    // Loads genetic map files required for recombination rate estimation
    genMap_ch = Channel.fromPath("${params.Relate.genMap_dir.replaceAll('/+$', '')}/*.genetic_map.txt")
    genMap_ch.ifEmpty { error "❌ No .genetic_map.txt files found in ${params.Relate.genMap_dir}" }

    // Extract chromosome information from genetic map files
    genMap_ch.map { file -> 
        def m = file.getName() =~ /^(?:(.*)[._])?chr([^\.]+)\.genetic_map\.txt$/
        if (!m) {
            error "Invalid filename format: ${file.name}"
        }else{
            tuple(m[0][2], file)
        }
    }.set { relate_genetic_map_files }
    // Process hap and sample files and combine with genetic maps
    relate_files
        .map { relate_sample, relate_hap -> 
            def m = relate_sample.getName() =~ /^(.*)_([^\.]+)\.sample\.gz$/
            if (!m) error "Invalid filename format: ${relate_sample.name}"
            tuple(m[0][2], relate_sample, relate_hap)
        }
        .set { relate_gMap_files }
    // Combine genetic maps with annotation files by chromosome
    relate_genetic_map_files
        .combine(relate_gMap_files, by:0)
        .set{genMap_relate_files}

    // Creates ancestral trees for each chromosome
    //make chunk
    makeChunks(genMap_relate_files,Relate_path).set{chr_chunks}
    chr_chunks.chunk_files
        .flatten()
        .map { file ->
            def match = file.toString() =~ /Relate_chr([^\/]+)\/(?:chunk_|parameters_c)(\d+)\./
            if (match.find()) {
                def chrom = match.group(1)
                def chunk_id = match.group(2)
                return tuple(chrom, chunk_id, file)
            } else {
                return null
            }
        }
        .filter { it != null }
        .groupTuple(by: [0, 1], size: 7)
        .set { grouped_chunk_files }
    
    painting = '0.025,1'
    batch_windows = 5
    seed = 1
    mutation_rate = params.mutation_rate
    initial_N = params.Relate.initial_N
    paintChunks(grouped_chunk_files,painting,Relate_path).set{painted_chunks}
    getWindows(painted_chunks).set{windows_number}
    
    windows_number.flatMap { chr, chunk, file -> 
        def matcher = file.getName() =~ /^window_num_([^\.]+)$/
        def window_num = matcher[0][1].toInteger()
        def num_batched_windows = window_num / batch_windows
        return (0..<num_batched_windows).collect{ x -> tuple(chr, chunk, x) } }
        .set { index_chunk_window }

    buildTopology(index_chunk_window,batch_windows,seed,painting,Relate_path).set{tree_topology}
    tree_topology
        .map { chr, chunk, window, chunk_file -> def matcher = chr =~ /^Relate_chr([^\.]+)$/
        return tuple(matcher[0][1], chunk, window, chunk_file)}
        .groupTuple(by: [0,1]).combine(grouped_chunk_files,by:[0,1])
        .map { chr, chunk, grouped_windows, chunk_file, chunk_paras ->
        tuple(
            chr,
            chunk,
            chunk_file[0],
            grouped_windows,
            chunk_paras // keep chunk_paras unflattened
        )
    }.set { tree_topology_grouped }
    findEquivalentBranches(tree_topology_grouped,Relate_path).set { eq_branches }
    eq_branches.flatMap{
        chr, params, chunk, all_windows -> 
        all_windows.collect{ x -> tuple(chr, chunk, x)}}
        .set{ eq_branches_chunk_windows }
    inferBranchLengths(eq_branches_chunk_windows,batch_windows,mutation_rate,initial_N,seed,Relate_path).set{inferred}
    inferred
        .map{dir, chunk, window -> def matcher = dir.getName() =~ /^Relate_chr([^\.]+)$/
        return tuple(matcher[0][1], chunk, dir, window) }
        .groupTuple(by: [0,1]).map{
            chr, chunk, dir, window ->
            tuple(dir[0],chunk)
        }
        .set{inferred_all_windows}
    combineArgs(inferred_all_windows,initial_N,Relate_path).set{args_dir}
    chr_chunks.parameters_props.map{ parameters, props->
        def matcher = (parameters.toString() =~ /Relate_chr(\d+)[\/.]/)
        tuple(matcher[0][1], parameters, props)
    }
    .set { chr_params_props }
    args_dir.map{
        dir, chunks -> def matcher = dir.getName() =~ /^Relate_chr([^\.]+)$/
        return tuple(matcher[0][1],chunks)}
        .groupTuple(by: [0]).combine(chr_params_props,by: 0)
        .map{
            chr, chunks, parameters, props  ->
            tuple(chr,parameters,props, chunks.flatten())
        }
        .set{args_all_chunks}
    relate_tree_files = finalizeRelate(args_all_chunks,Relate_path)
    chr = relate_tree_files.chr_value
    relate_file = relate_tree_files.relate_files

}

process split_ancestral_chr{
    tag "split ancestral genome"
    input:
    path(ancestor_fast)

    output:
    path("*.fa")

    script:
    """
    # Function to process FASTA content
    process_fasta() {
        local current_chr=""
        local current_file=""

        while IFS= read -r line; do
            if [[ \$line == \\>* ]]; then
                # Header line - extract chromosome name
                chr_name=\$(echo "\$line" | sed 's/^>//' | awk '{print \$1}')
            
                # Close previous file if open
                if [ -n "\$current_file" ]; then
                    exec 3>&-
                fi
            
                # Open new file for this chromosome
                current_chr="\$chr_name"
                current_file="\${chr_name}.fa"
            
                echo "Processing: \$chr_name -> \${chr_name}.fa"
            
                # Open file descriptor for writing
                exec 3>"\$current_file"
                echo "\$line" >&3
            else
                # Sequence line - write to current file
                if [ -n "\$current_file" ]; then
                    echo "\$line" >&3
                fi
            fi
        done
    
        # Close last file
        if [ -n "\$current_file" ]; then
            exec 3>&-
        fi
    }
    if [[ ${ancestor_fast} == *.gz ]]; then
        echo "Processing gzipped FASTA file..."
        zcat ${ancestor_fast} | process_fasta
    else
        echo "Processing uncompressed FASTA file..."
        process_fasta < ${ancestor_fast}
    fi

    """
}

process prepare_relate_format {
    tag "Preparing Relate format: ${vcf_file.baseName}"
    stageInMode = 'symlink'
    input:
    tuple(val(chr),path(selected_poplabel), path(vcf_file), val(basename), path(ancestor))  // Input VCF file and population labels
    val(Relate_path)  // Path to Relate software

    output:
    tuple(path("${basename}.sample.gz"), path("${basename}.haps.gz"))  // Output files in Relate format

    script:
    """
    # check Relate path
    if [ ! -d "${params.Relate.Relate_path}" ]; then
        echo "ERROR: Relate_path not found: ${params.Relate.Relate_path}"
        exit 1
    fi

    # check Relate executable
    if [ ! -f "${params.Relate.Relate_path}/bin/RelateFileFormats" ]; then
        echo "ERROR: Relate executable not found at: ${params.Relate.Relate_path}/bin/RelateFileFormats"
        ls ${params.Relate.Relate_path} || echo "Directory is empty or inaccessible"
        exit 1
    fi

    name_vcf=\$(basename "${vcf_file}" .vcf.gz)
    tail -n +2 ${selected_poplabel} | cut -f1 > known_samples.txt
    zgrep '^#CHROM' "${vcf_file}" | cut -f10- | tr '\\t' '\\n' > vcf_samples.txt
    comm -23 <(sort vcf_samples.txt) <(sort known_samples.txt) > missing.txt

    # Convert VCF to Relate format
    ${Relate_path}/bin/RelateFileFormats \\
        --mode ConvertFromVcf \\
        --haps ${basename}.haps \\
        --sample ${basename}.sample \\
        -i \${name_vcf}
    ${Relate_path}/bin/RelateFileFormats \\
        --mode RemoveNonBiallelicSNPs \\
        --haps ${basename}.haps \\
        -o ${basename}_biall

    ${Relate_path}/bin/RelateFileFormats \\
        --mode FlipHapsUsingAncestor \\
        --haps ${basename}_biall.haps \\
        --sample ${basename}.sample \\
        --ancestor ${ancestor} \\
        -o ${basename}_ancest
    rm ${basename}_biall.haps

    ${Relate_path}/bin/RelateFileFormats \\
        --mode RemoveSamples \\
        --haps ${basename}_ancest.haps \\
        --sample ${basename}.sample \\
        -i missing.txt \\
        -o ${basename}_removed

    rm ${basename}_ancest.haps
    mv ${basename}_removed.haps ${basename}.haps
    mv ${basename}_removed.sample ${basename}.sample
    gzip ${basename}.haps
    gzip ${basename}.sample
    """
}

process makeChunks {
    tag "make chunk"
    stageInMode = 'symlink'
    input:
    tuple(val(chr),path(genMap),path(relate_sample), path(relate_hap))
    val(Relate_path)

    output:
    val(chr), emit: chr
    path("Relate_chr${chr}/{chunk_*,parameters_c*}"), emit: chunk_files
    path("Relate_chr${chr}/{parameters.bin,props.bin}"),emit: parameters_props

    script:
    """
    ${Relate_path}/bin/Relate \\
        --mode "MakeChunks" \\
        --haps ${relate_hap} \\
        --sample ${relate_sample} \\
        --map ${genMap} \\
        --memory 5 \\
        --output Relate_chr${chr}
    """
}

process paintChunks {
    tag "paint Relate_chr${chr}_${chunk_index}"
    stageInMode = 'symlink'
    input:
    tuple(val(chr),val(chunk_index),path(chunk_files,stageAs: "Relate_chr/"))
    val(painting)
    val(Relate_path)

    output:
    tuple(path("Relate_chr${chr}/"),val(chunk_index),path("Relate_chr${chr}/parameters_c${chunk_index}.bin"))

    script:
    """
    mkdir -p Relate_chr${chr}/
    mv Relate_chr/* Relate_chr${chr}/
    ${Relate_path}/bin/Relate \\
    --mode "Paint" \\
    --chunk_index ${chunk_index} \\
    --painting ${painting} \\
    --output Relate_chr${chr}

    # generate pain/*.bin
    """
}

process getWindows {
    tag "get windows from ${chunk_dir}_${chunk_index}"
    stageInMode = 'symlink'
    label "basic_Ne"
    input:
        tuple(path(chunk_dir),val(chunk_index),path(parameters_bin))
    output:
        tuple(path(chunk_dir),val(chunk_index),path("window_num_*"))

    script:
    """
    ACTIVE_ENV=\$(conda info | grep 'active environment' | awk '{print \$NF}')
    if [ "\$ACTIVE_ENV" = "None" ] || [ "\$ACTIVE_ENV" != "basic_Ne" ]; then
        source \$(conda info --base)/etc/profile.d/conda.sh
        conda activate basic_Ne
    fi
    num_windows=(\$(python -c "
    import numpy as np
    a = np.fromfile('${parameters_bin}', dtype=np.int32)
    print(a[0], a[1], a[2])
    "))
    num_windows=\${num_windows[2]}
    echo \$num_windows > window_num_\${num_windows}
    """
}
process buildTopology {
    tag "Build Topology: chunk-${chunk_index} windows-${batch_id}"
    stageInMode = 'symlink'
    input:
    tuple(path(folder),val(chunk_index),val(batch_id))
    val(batch_windows)
    val(seed)
    val(painting)
    val(Relate_path)

    output:
    tuple(val("${folder}"),val("${chunk_index}"), val("${batch_id}"), path("${folder}/chunk_${chunk_index}/"))

    script:
    """
    start_ind=\$(( ${batch_id} * ${batch_windows} ))
    end_ind=\$(( (${batch_id} + 1) * ${batch_windows} - 1 ))
    ${Relate_path}/bin/Relate \\
    --mode "BuildTopology" \\
    --first_section \$start_ind \\
    --last_section \$end_ind \\
    --chunk_index ${chunk_index} \\
    --painting ${painting} \\
    -o ${folder}

    # generate .chunk_/Relate_chr*_chunk.anc and .mut
    """
}

process findEquivalentBranches {
    tag 'Find Equivalent Branches'
    stageInMode = 'copy'
    input:
    tuple(val(chr),val(chunk_index),path(chunk_files),val(all_windows),path(chunk_paras))
    val(Relate_path)
    output:
    tuple(path("Relate_chr${chr}"),path("Relate_chr${chr}/*.r"),val(chunk_index),val(all_windows))
    
    script:
    """
    mkdir -p Relate_chr${chr}/
    mv ${chunk_files} Relate_chr${chr}/
    mv ${chunk_paras} Relate_chr${chr}/
    # Remove paint directory contents if it exists
    if [ -d "Relate_chr${chr}/paint" ]; then
        rm -rf Relate_chr${chr}/paint/*
    fi

    ${Relate_path}/bin/Relate \\
    --mode "FindEquivalentBranches" \\
    --chunk_index ${chunk_index} \\
    -o Relate_chr${chr}
    """
}

process inferBranchLengths {
    tag 'Infer Branch Lengths'
    input:
    tuple(path(topo_infered_dir),val(chunk_index),val(window_index))
    val(batch_windows)
    val(mutation_rate)  // Mutation rate parameter
    val(initial_N)  // Initial population size
    val(seed)
    val(Relate_path)

    output:
    tuple(path(topo_infered_dir),val(chunk_index),val(window_index))

    script:
    """
    start_ind=\$(( ${window_index} * ${batch_windows} ))
    end_ind=\$(( (${window_index} + 1) * ${batch_windows} - 1 ))
    ${Relate_path}/bin/Relate \\
    --mode "InferBranchLengths" \\
    -m ${mutation_rate} \\
    -N ${initial_N} \\
    --first_section \$start_ind \\
    --last_section \$end_ind \\
    --chunk_index ${chunk_index} \\
    --seed ${seed} \\
    -o ${topo_infered_dir}
    # delete .hap .r .rpos .state
    """
}

process combineArgs {
    tag 'Combine Args'
    stageInMode = 'symlink'
    input:
    tuple(path(args_dir), val(chunk_index))
    val(initial_N)  // Initial population size
    val(Relate_path)
    output:
    tuple(path(args_dir),path("${args_dir}/chunk_${chunk_index}"))
    script:
    """
    ${Relate_path}/bin/Relate \\
    --mode "CombineSections" \\
	-N ${initial_N} \\
	--chunk_index ${chunk_index} \\
	-o ${args_dir}

    """
}

process finalizeRelate {
    tag 'Finalize'
    publishDir "${params.output_dir}/Relate", mode: 'copy'
    stageInMode = 'copy'
    input:
    tuple(val(chr),path(parameter),path(probs),path(files))
    val(Relate_path)
    output:
    val(chr), emit: chr_value
    path("Relate_chr*"), emit: relate_files  // Output tree and mutation files

    script:
    """
    mkdir -p Relate_chr${chr}/
    mv ${parameter} ${probs} ${files} Relate_chr${chr}/
    ${Relate_path}/bin/Relate \\
    --mode "Finalize" \\
    -o Relate_chr${chr}
    """
}

process RemoveTrees {
    tag "Remove Trees With Few Mutations ${chr}"
    input:
    val Relate_path
    val output_prefix
    val threshold
    tuple val(chr), path(anc_file), path(mut_file)

    output:
    tuple val(chr), path("*.anc.gz"), path("*.mut.gz"), emit: anc_mut
    tuple val(chr), path("*.dist"), emit: dist
    
    script:
    """
    #delete all trees that have fewer than threshold mutations
    ${Relate_path}/bin/RelateExtract \\
    --mode RemoveTreesWithFewMutations \\
    --threshold ${threshold} \\
    --anc ${anc_file} \\
    --mut ${mut_file} \\
    -o ${output_prefix}_chr${chr}

    gzip ${output_prefix}_chr${chr}.anc
    gzip ${output_prefix}_chr${chr}.mut
    """
}


process SampleBranchCallCoal{
    input:
    val(Relate_path)
    val(output_prefix)
    val(iteration)
    val(chunk_num)
    val(mutation)
    path(anc_mut)
    path(dist)
    path(coal)
    path(chr_list)

    output:
    path("*.anc.gz"), emit: anc
    path("*.mut.gz"), emit: mut
    path("*.coal"), emit: coal

    script:
    seed_arg = params.Relate.seed ? "--seed ${params.Relate.seed}" : ""
    bins_arg = params.Relate.bins ? "--bins ${params.Relate.bins}" : ""
    years_per_gen = params.Relate.years_per_gen ?: 5
    """
    # export PARALLEL_HOME=/tmp/parallel_\$\$
    # export TMPDIR=/tmp
    # mkdir -p \$PARALLEL_HOME
    for i in \$(seq 1 ${iteration}); do
        # Step 1: Divide AncMut files for each chromosome in parallel
        cat ${chr_list} | xargs -P ${task.cpus} -I {} sh -c '
            chr=\$1
            ${Relate_path}/bin/RelateExtract --mode DivideAncMut --threads ${chunk_num} --anc ${output_prefix}_chr\${chr}.anc.gz --mut ${output_prefix}_chr\${chr}.mut.gz -o ${output_prefix}_chr\${chr}_tmp
        ' _ {}

        # Step 2: Process all chromosome-chunk combinations
        # First, get all chromosomes and their chunk counts
        cat ${chr_list} | while read chr; do
            count=\$(ls ${output_prefix}_chr\${chr}_tmp_chr*.mut.gz 2>/dev/null | wc -l)
            
            seq 0 \$((\$count - 1)) | xargs -P ${task.cpus} -I {} sh -c '
                chr=\$1
                benchmark=\$2
                ${Relate_path}/bin/RelateCoalescentRate --mode SampleBranchLengths --coal ${output_prefix}.coal --dist ${output_prefix}_chr\${chr}.dist --num_samples 1 ${seed_arg} -m ${mutation} -i ${output_prefix}_chr\${chr}_tmp_chr\${benchmark} -o ${output_prefix}_chr\${chr}_tmp_chr\${benchmark} && rm ${output_prefix}_chr\${chr}_tmp_chr\${benchmark}.anc.gz && rm ${output_prefix}_chr\${chr}_tmp_chr\${benchmark}.mut.gz
            ' _ "\$chr" {}
        done
        # Step 3: Combine AncMut files and cleanup
        cat ${chr_list} | xargs -P ${task.cpus} -I {} sh -c '
            chr=\$1
            ${Relate_path}/bin/RelateExtract --mode CombineAncMut -o ${output_prefix}_chr\${chr}_tmp
            rm ${output_prefix}_chr\${chr}.anc.gz
            rm ${output_prefix}_chr\${chr}.mut.gz
            mv ${output_prefix}_chr\${chr}_tmp.anc.gz ${output_prefix}_chr\${chr}.anc.gz
            mv ${output_prefix}_chr\${chr}_tmp.mut.gz ${output_prefix}_chr\${chr}.mut.gz
        ' _ {}
        
        ${Relate_path}/bin/RelateCoalescentRate \\
            --mode CoalRateForTree \\
            --years_per_gen ${years_per_gen} \\
            ${bins_arg} \\
            --chr ${chr_list} \\
            -i ${output_prefix} \\
            --dist ${output_prefix} \\
            -o ${output_prefix}_tmp
        rm ${output_prefix}.coal
        mv ${output_prefix}_tmp.coal ${output_prefix}.coal

    done

    cat ${chr_list} | while read chr; do
        mv ${output_prefix}_chr\${chr}.anc.gz ${output_prefix}_final_chr\${chr}.anc.gz
        mv ${output_prefix}_chr\${chr}.mut.gz ${output_prefix}_final_chr\${chr}.mut.gz
    done
    mv ${output_prefix}.coal ${output_prefix}_final.coal
    """

}



/*
 * Process: CoalRateForTree
 * Purpose: Estimate genome-wide coalescence rate across all chromosomes.
 *          Produces a .coal file used as the starting point for the
 *          SampleBranchCallCoal iteration loop.
 * Input:  chr list, anc/mut files, dist files
 * Output: *.chr (chromosome list), *.coal (coalescence rate file)
 */
process CoalRateForTree {
    publishDir "${params.output_dir}/Relate"
    input:
    val(Relate_path)
    val(output_prefix)
    val(chr_list)
    path(anc_mut)
    path(dist)

    output:
    path("*.chr"), emit: chr_list
    path("*.coal"), emit: coal

    script:
    bins_arg = params.Relate.bins ? "--bins ${params.Relate.bins}" : ""
    years_per_gen = params.Relate.years_per_gen ?: 5
    """
    echo ${chr_list} | tr -d '[] ' | tr ',' '\\n' | sed 's/^chr//' > ${output_prefix}.chr
    ${Relate_path}/bin/RelateCoalescentRate \\
        --mode CoalRateForTree \\
        --years_per_gen ${years_per_gen} \\
        ${bins_arg} \\
        --chr ${output_prefix}.chr \\
        --dist ${output_prefix} \\
        -i ${output_prefix} \\
        -o ${output_prefix}
    """
}

/*
 * Process: CoalRatePairwise
 * Purpose: Compute pairwise coalescence rates per chromosome section.
 *          These per-chromosome results are later aggregated by
 *          SummarizeCoalRate into a genome-wide Ne estimate.
 * Input:  chr, anc/mut/dist/coal files
 * Output: *.pairwise_chr<N> files
 */
process CoalRatePairwise {
    input:
    val(Relate_path)
    val(output_prefix)
    path(chr_list)
    tuple(val(chr),path(anc_file), path(mut_file), path(dist_file), path(coal_file))

    output:
    path("*.pairwise_chr*")

    script:
    bins_arg = params.Relate.bins ? "--bins ${params.Relate.bins}" : ""
    years_per_gen = params.Relate.years_per_gen ?: 5
    """
    ${Relate_path}/bin/RelateCoalescentRate \\
        --mode CoalescentRateForSection \\
        --years_per_gen ${years_per_gen} \\
        ${bins_arg} \\
        --dist ${dist_file} \\
        --chr ${chr_list} \\
        -i ${output_prefix}_chr${chr} \\
        -o ${output_prefix}.pairwise_chr${chr}
    """
}

process SummarizeCoalRate {
    publishDir "${params.output_dir}/Relate", mode: 'copy'
    input:
    val(Relate_path)
    val(output_prefix)
    path(chr_file)
    path(pairwise_files)
    path(poplabels_file)

    output:
    path("*.pairwise.coal")
    
    script:
    bins_arg = params.Relate.bins ? "--bins ${params.Relate.bins}" : ""
    years_per_gen = params.Relate.years_per_gen ?: 5
    """
    # Summarize coalescent rate for genome
    ${Relate_path}/bin/RelateCoalescentRate \\
        --mode SummarizeCoalescentRateForGenome \\
        --chr ${chr_file} \\
        -o ${output_prefix}.pairwise
    
    # Finalize population size
    ${Relate_path}/bin/RelateCoalescentRate \\
        --mode FinalizePopulationSize \\
        -o ${output_prefix}.pairwise \\
        --poplabels ${poplabels_file}
    
    """
}

/*
 * Process: DivideAncMutFinal
 * Purpose: Same as DivideAncMut but applied to the final iterated anc/mut
 *          files before ReEstimateBranchLengths.  Also emits a .param file
 *          required by CombineAncMutFinal.
 * Input:  chr, final anc/mut files
 * Output: per-chunk *_tmp_chr* files + *.param metadata
 */
process DivideAncMutFinal{
    input:
    val(Relate_path)
    val(output_prefix)
    val(chunk_num)
    tuple(val(chr),path(anc_file), path(mut_file))

    output:
    tuple val(chr), path("*_tmp_chr*"), emit: anc_mut
    tuple val(chr), path("*.param"),    emit: params_file

    script:
    """
    ${Relate_path}/bin/RelateExtract \\
        --mode DivideAncMut \\
        --threads ${chunk_num} \\
        --anc ${anc_file} \\
        --mut ${mut_file} \\
        -o ${output_prefix}_chr${chr}_tmp
    """

}

/*
 * Process: ReEstimateBranchLengths
 * Purpose: Re-estimate branch lengths per chunk using the final coalescence
 *          rates.  Analogous to SampleBranchLength but runs after the
 *          iterative loop has converged.
 * Input:  chr, chunk index, anc/mut chunks, dist file, coal file
 * Output: updated *_tmp_chr* chunk files
 */
process ReEstimateBranchLengths {
    input:
    val(Relate_path)
    val(output_prefix)
    val(mutation)
    tuple(val(chr),val(chunk), path(anc_mut), path(dist_file), path(coal_file))

    output:
    tuple(val(chr), path("*_tmp_chr*"))

    script:
    seed_arg = params.Relate.seed ? "--seed ${params.Relate.seed}" : ""
    """
    ${Relate_path}/bin/RelateCoalescentRate \\
        --mode ReEstimateBranchLengths \\
        --coal ${coal_file} \\
        --dist ${dist_file} \\
        -m ${mutation} \\
        ${seed_arg} \\
        -i ${output_prefix}_chr${chr}_tmp_chr${chunk} \\
        -o ${output_prefix}_chr${chr}_tmp_chr${chunk}
    """
}

/*
 * Process: CombineAncMutFinal
 * Purpose: Merge re-estimated per-chunk files into final per-chromosome
 *          anc/mut files.  These are the definitive Relate tree files
 *          published to the output directory.
 * Input:  chr, chunk files, .param metadata
 * Output: final *_chr*.anc.gz and *_chr*.mut.gz
 */
process CombineAncMutFinal {
    publishDir "${params.output_dir}/Relate", mode: 'copy'
    input:
    val(Relate_path)
    val(output_prefix)
    tuple(val(chr),path(files),path(params_file))

    output:
    tuple(path("*_chr*.anc.gz"), path("*_chr*.mut.gz"))

    script:
    """
    ${Relate_path}/bin/RelateExtract \\
        --mode CombineAncMut \\
        -o ${output_prefix}_chr${chr}_tmp
    mv ${output_prefix}_chr${chr}_tmp.anc.gz ${output_prefix}_chr${chr}.anc.gz
    mv ${output_prefix}_chr${chr}_tmp.mut.gz ${output_prefix}_chr${chr}.mut.gz
    """
}

process summarize_Relate {
    label "basic_Ne"
    tag "Summarizing Relate results"
    publishDir "${params.output_dir}", mode: 'copy'
    input:
    path(relate_output_files)  // Input files containing Ne estimates
    
    output:
    path("Relate_results.json")  // Output JSON file with summary statistics
    
    script:
    """
    ACTIVE_ENV=\$(conda info | grep 'active environment' | awk '{print \$NF}')
    if [ "\$ACTIVE_ENV" = "None" ] || [ "\$ACTIVE_ENV" != "basic_Ne" ]; then
        source \$(conda info --base)/etc/profile.d/conda.sh
        conda activate basic_Ne
    fi

    python ${projectDir}/scripts/summarize_results.py \\
        -m "Relate" \\
        -i ${relate_output_files} \\
        -o Relate_results.json
    """
}




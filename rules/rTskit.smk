configfile: "./config/tskit_config.yaml"

rule run_prepare_ts:
    input:
        vcf = f'{workDir}/{vcf_f}/{{file_number}}.vcf.gz'
        ancestal = config["ancestral_genome"]
        popLabel = f'{workDir}/{vcf_f}/{{file_number}}.popLabel'
    output:
        # Expanding OUTPUT_data for each sampling number
        expand(f'{workDir}/{vcf_f}/tskit_results/{{file_number}}.samples',file_number='{file_number}')
    conda: "HLab_tsinfer"
    script:
        "./script/ts_infer.py"
        

rule run_tsinfer:
    input:
        f'{workDir}/{vcf_f}/tskit_results/{{sample_file}}.samples'
    output:
        f'{workDir}/{vcf_f}/tskit_results/{{sample_file}}.trees'
    conda: "HLab_tsinfer"
    run:
        ts = tsinfer.infer(input)
        print(
            "Inferred tree sequence: {} trees over {} Mb ({} edges)".format(
        ts.num_trees, ts.sequence_length / 1e6, ts.num_edges
                )
                    )
        redated_ts = tsdate.date(ts, mutation_rate=args.mutation)
        redated_ts.dump(output)


rule summary_tsinfer:
    input:
        ts_folder = f'{workDir}/{vcf_f}/tskit_results/'
        timepoint = config["timepoint"]
        output = f'{workDir}/{vcf_f}/tskit_results.json'
    output:
        f'{workDir}/{vcf_f}/tskit_results.json'
    script:
        "./script/ts_summary_Ne.py"
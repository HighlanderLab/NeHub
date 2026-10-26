#!/usr/bin/env python
import os
import pysam
import phlash
import json
import numpy as np
import pickle
import argparse
import subprocess

def get_last_vcf_record(vcf_path):
    try:
        result = subprocess.run(
            f"zgrep -v '^#' {vcf_path} | tail -n 1",
            shell=True,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True
        )
        line = result.stdout.strip()
        if not line:
            return None
        fields = line.split('\t')
        return int(fields[1])
    except subprocess.CalledProcessError:
        return None
def run_phlash(input_vcf_phlash, input_sampleid, output_file, mutation_rate,chr_test):
    """
    Run the Phlash analysis using the specified VCF file, population labels, and mutation rate.

    Parameters:
        input_vcf_phalsh (str): Path to the VCF file.
        input_sampleid (list of str): List of sample.
        output_dir (str): Directory to save the results.
        mutation_rate (float): Mutation rate for the Phlash fit.

    Returns:
        None
    """


    # Get sample IDs for the population
    samples_list = input_sampleid
    print(samples_list)
    # Analyze each chromosome for the population
    chroms_1kg = []
    for input_vcf in input_vcf_phlash:
        # Open the VCF file
        vcf = pysam.VariantFile(input_vcf)
        chrom = list(vcf.header.contigs)
        chrom = chrom[0]
        length = vcf.header.contigs[chrom].length
        if length == None:
            last_record = get_last_vcf_record(input_vcf)
            region = f"{chrom}:1-{last_record+1}"
        else:
            region = f"{chrom}:1-{length}"

        print(region)
        # Add the contig analysis to the list
        chroms_1kg.append(
            phlash.contig(input_vcf, samples=samples_list, region=region)
            )
    print("start fit model")
    # Fit the results using Phlash
    if chr_test==0:
        results = phlash.fit(chroms_1kg, mutation_rate=mutation_rate)
    else:
        test_data = chroms_1kg[0]
        train_data = chroms_1kg[1:]
        results = phlash.fit(data=train_data, test_data=test_data, mutation_rate=mutation_rate)
    # Save results to a pickle file
    with open(output_file+'.pkl', "wb") as f:
        pickle.dump(results, f)
    
    with open(output_file+'.json', "w") as f:
        times = np.array([dm.eta.t[1:] for dm in results])
        # choose a grid of points at which to evaluate the size history functions
        T = np.geomspace(times.min(), times.max(), 100)
        Nes = np.array([dm.eta(T, Ne=True) for dm in results])
        T_Nes = {}
        for index in range(len(T)):
            Ne = Nes[:,index]
            T_Nes[str(float(T[index]))] = [
                float(np.mean(Ne)), 
                float(np.median(Ne)), 
                float(np.percentile(Ne, 2.5)), 
                float(np.percentile(Ne, 97.5)), 
                float(np.percentile(Ne, 12.5)), 
                float(np.percentile(Ne, 87.5))
            ]
        json.dump(T_Nes, f, indent=4)
    print("Phlash analysis complete. Results saved.")

# Example of how to use the function in a Snakemake workflow
if __name__ == "__main__":
    # Create an ArgumentParser object
    parser = argparse.ArgumentParser(description="Process some integers.")

# Add arguments
    parser.add_argument('--vcfs', nargs="+", required=True, help='Path to the VCF file')
    parser.add_argument('--sampleId', nargs="+", required=True, help='sample Id')
    parser.add_argument('--mutation_rate', type=float, required=True,help='mutation rate')
    parser.add_argument('--chr_test', type=int, default=0, help='if set the test data, if yes, the first chr is test data')
    parser.add_argument('--out_file',type=str, required=True, help='output file')
    # Parse arguments
    args = parser.parse_args()
    print(args.sampleId)
    run_phlash(
        input_vcf_phlash=args.vcfs,
        input_sampleid=args.sampleId,
        output_file=args.out_file,
        mutation_rate=args.mutation_rate,
        chr_test=args.chr_test
        )

import pandas as pd
import os
import pysam
import phlash
import pickle
import argparse
def run_phlash(input_vcf_phalsh, input_pop_labels, output_dir, mutation_rate):
    """
    Run the Phlash analysis using the specified VCF file, population labels, and mutation rate.

    Parameters:
        input_vcf_phalsh (str): Path to the VCF file.
        input_pop_labels (list of str): List of file paths for population labels.
        output_dir (str): Directory to save the results.
        mutation_rate (float): Mutation rate for the Phlash fit.

    Returns:
        None
    """

    # Change to the working directory
    os.chdir(output_dir)

    # Open the VCF file
    vcf = pysam.VariantFile(input_vcf_phalsh)

    # Loop through input population labels
    
    df = pd.read_csv(input_pop_labels, sep="\t")  # Read population data
    pops = df['population'].unique()         # Get unique populations

    for pop in pops:
        # Get sample IDs for the population
        samples_list = df.loc[df['population'] == pop, 'sample_id'].to_list()

        # Analyze each chromosome for the population
        chroms_1kg = []
        for chrom in vcf.header.contigs:
            length = vcf.header.contigs[chrom].length
            region = f"{chrom}:1-{length}"

            # Add the contig analysis to the list
            chroms_1kg.append(
                    phlash.contig(input_vcf_phalsh, samples=samples_list, region=region)
                )

        # Fit the results using Phlash
        results = phlash.fit(chroms_1kg, mutation_rate=mutation_rate)
        # Save results to a pickle file
        output_file = os.path.join(output_dir, f"{pop}.pkl")
        with open(output_file, "wb") as f:
            pickle.dump(results, f)

    print("Phlash analysis complete. Results saved.")

# Example of how to use the function in a Snakemake workflow
if __name__ == "__main__":
    # Create an ArgumentParser object
    parser = argparse.ArgumentParser(description="Process some integers.")

# Add arguments
    parser.add_argument('--vcf', type=str, required=True, help='Path to the VCF file')
    parser.add_argument('--poplabels', type=str, required=True, help='Path to the ancestral genome fa file')
    parser.add_argument('--mutation_rate', type=float, required=True,help='mutation rate')
    parser.add_argument('--out',type=str, required=True, help='output folder')
    # Parse arguments
    args = parser.parse_args()
    run_phlash(
        input_vcf_phalsh=args.vcf,
        input_pop_labels=args.poplabels,
        output_dir=args.out,
        mutation_rate=args.mutation_rate
        )

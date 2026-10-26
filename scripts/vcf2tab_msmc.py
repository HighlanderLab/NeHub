#!/usr/bin/env python
"""
Convert one or more VCF files to the MSMC2 tab-delimited input format.

For each population×group in the population file the script:
  1. Randomly samples up to sampleSize diploid individuals.
  2. For each input VCF extracts phased alleles for those individuals.
  3. Writes one <pop>_<group>_<chr>_<seed>.tab file per chromosome.
  4. Also writes a <pop>_<group>_<seed>.csv with the sampled metadata.

MSMC tab format (per line):
    CHROM  POS  INTERVAL  ALLELES…
where ALLELES is the concatenation of all haplotype bases (e.g. "AACGTG").

Usage:
    python vcf2tab_msmc.py -i chr1.vcf.gz chr2.vcf.gz \
        -p popLabel.txt -s 4
"""

import os
import sys
import argparse
import pandas as pd
import cyvcf2
import time


def readvcf(vcf, pop_list):
    """Process VCF records for current sample group"""
    dic_output = {}
    prev_pos = 0  # Track previous position for interval calculation
    
    for record in vcf:
        if prev_pos==0:
            i_index = [vcf.samples.index(s) for s in pop_list]
        chrom = record.CHROM
        position = record.POS
        
        # Calculate interval from previous position
        interval = position - prev_pos
        prev_pos = position

        # Process genotypes for all samples in current group
        alleles_str = ""
        for gt in record.gt_bases[i_index]:
            if (gt is None) or (gt == "./."):
                # Handle missing genotypes
                alleles_str += "NN"
                continue
                
            # Split phased/unphased genotypes and join alleles
            alleles = gt.replace('|', '/').split('/')
            if len(alleles) != 2:
                # Handle non-diploid or malformed genotypes
                alleles_str += "NN"
                continue
                
            alleles_str += ''.join(alleles)

        # Store formatted line
        line = f"{chrom}\t{position}\t{interval}\t{alleles_str}\n"
        dic_output.setdefault(chrom, []).append(line)
    
    return dic_output

def vcf2tab(vcf_file, pop_name, pop_list,seed):

    print(f"Processing group {pop_name} with {len(pop_list)} samples ")
    print(",".join(pop_list))
    try:
        # Create new VCF reader for each group with subset of samples
        vcf = cyvcf2.VCF(vcf_file)
        dic_chr=readvcf(vcf, pop_list)
        # Read population metadata
        for chrom, lines in dic_chr.items():
            outfile = f"{pop_name}_{chrom}_{seed}.tab"
            with open(outfile, "w") as f:
                f.writelines(lines)
    except Exception as e:
        sys.exit(f"Error processing group {pop_name}: {str(e)}")


def read_population(file_path):
    """Read population metadata file"""
    if not os.path.exists(file_path):
        sys.exit(f"Error: File not found: {file_path}")
    return pd.read_csv(file_path, sep="\t")

def main():
    parser = argparse.ArgumentParser(description="Process VCF files by population groups")
    parser.add_argument("-i", "--input", required=True, nargs='+', help="Input VCF file")
    parser.add_argument("-p", "--pop_file", required=True, help="Population file to process")
    parser.add_argument("-s", "--sampleSize", required=True, help="Sampling size, default=4 (4 individuals) MSMC has limitation of at most 4 individuals per group")
    parser.add_argument("-seed", "--seed", required=True, help="Seed for random number generator")
    args = parser.parse_args()
    vcf_files = args.input
    pop_file = args.pop_file
    sample_size = int(args.sampleSize)
    seed = int(args.seed)
    population=read_population(pop_file)
    if 'group' not in population.columns:
        population["group"]=population["population"]
    
    if 'sex' not in population.columns:
        population['sex'] = "NA"

    # Perform stratified sampling on population data
    groups = [group for name, group in population.groupby(["population","group"])]
    for group in groups:
        df_group=group.sample(sample_size,random_state=seed)
        df_group.to_csv(f"{df_group['population'].values[0]}_{df_group['group'].values[0]}_{seed}.csv", index=False)
        for  vcf_file in vcf_files:
            # Process each VCF file for the current sample group
            if not os.path.exists(vcf_file):
                sys.exit(f"Error: VCF file not found: {vcf_file}")
            else:
                vcf2tab(vcf_file, f"{df_group['population'].values[0]}_{df_group['group'].values[0]}", df_group["sample_id"].tolist(),seed)
            # Call the function to process the VCF file
        

if __name__ == "__main__":
    main()
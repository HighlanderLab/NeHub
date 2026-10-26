#!/usr/bin/env python
"""
Stratified sampler for the Ne-pipeline population label file.

Reads a TSV population file (columns: sample_id, population[, group, sex])
and writes a new TSV containing a stratified sample of individuals.

Sampling behaviour (--sampling / -s):
  0 < s ≤ 1   – fraction of each population×group to keep (e.g. 0.5 = 50 %)
  s > 1        – fixed count per population×group (rounded, minimum 1)
  s == 1       – keep all (no sampling)

Usage:
    python sampler.py -i popLabel.txt -o selected_poplabel1.txt -s 1
"""

import argparse
import sys
import os
import pandas as pd
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description="Population File Sampler")
    parser.add_argument("-i", "--input", required=True, help="Input population file")
    parser.add_argument("-o", "--output", required=True, help="Output file path")
    parser.add_argument("-s", "--sampling", type=float, default=1,
                    help="Sampling rate (0.0-1.0), default=1.0 (no sampling), if the rate is greater than 1, it will be considered as the number of samples to be taken from each group")
    parser.add_argument("-seed", "--seed", type=int, default=42,
                    help="Seed for random number generator")
    return parser.parse_args()

def read_population(file_path):
    # Read population data
    if not os.path.exists(file_path):
        sys.exit(f"Error: File not found: {file_path}")
    population = pd.read_csv(file_path, sep="\t")
    population.rename(columns={population.columns[0]: 'sample_id',population.columns[1]: 'population'}, inplace=True)
    return population

def stratified_sampling(data, rate, seed):    
    
    sampled = pd.DataFrame()
    groups = [group for name, group in data.groupby(["population","group"])]
    for group in groups:
        if rate >1:
            n=rate
        else:
            n = max(1, int(len(group) * rate))
        sampled = pd.concat([sampled,group.sample(int(n),random_state=seed)],axis=0)
        
    return sampled


def main():
    args = parse_args()
    
    # Validate sampling rate
    if not 0 <= args.sampling:
        sys.exit(f"Invalid sampling rate: {args.sampling}. Must be >0")
    else:
        population = read_population(args.input)
        if 'group' not in population.columns:
            population["group"]=population["population"]
    
        if 'sex' not in population.columns:
            population['sex'] = "NA"

        if args.sampling == 1.0:
            population[['sample_id','population','group','sex']].to_csv(args.output, sep="\t", index=False)
        else:
            sampled = stratified_sampling(population, args.sampling,args.seed)
            sampled[['sample_id','population','group','sex']].to_csv(args.output, sep="\t", index=False)
            
        
if __name__ == "__main__":
    main()
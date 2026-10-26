#!/usr/bin/env python
"""
Summarize Ne-estimation results from all supported methods into a JSON file.

Usage:
    python summarize_results.py -m <METHOD> -i <file1> [file2 ...] -o <output.json>

Supported methods and their expected input:
    GONE     – Output_Ne_* files (two-column TSV: generation, Ne)
    MSMC     – *.final.txt files (MSMC2 output with lambda column)
    SMCpp    – *.csv files produced by 'smc++ plot'
    stairway – *.final.summary files from Stairway Plot 2
    Phlash   – *.json files from phlash_run.py
    PSMC     – *.psmc files (requires -g population-label file)
    Relate   – *.pairwise.coal files from SummarizeCoalRate
"""

import json
import argparse
import numpy as np
import os
import re
import pandas as pd
import math
import sys

def summarize_GONE_results(input_files, output):
    combined = {}
    for input_file in input_files:
        with open(input_file, 'r') as f:
            _ = f.readline()  # Skip annotation
            _ = f.readline()  # Skip header
            dic = {}
            for line in f:
                if not line.strip():
                    continue
                print(line)
                col1, col2 = line.strip().split('\t')
                dic[int(col1)] = float(col2)
            combined[input_file] = dic
    with open(output, 'w') as out_file:
        json.dump(combined, out_file, indent=4)
def summarize_MSMC_results(msmc_results, output, mutation_rate):
    combined = {}
    for i in msmc_results:
        with open(i, "r") as f:
            content = f.readlines()
        Ne = []
        for j in content[1:]:
            li_j = j.strip('\n').split('\t')
            t = float(li_j[1])
            la = float(li_j[3])
            Ne.append([round(t / mutation_rate), (1.0 / la) / (2.0 * mutation_rate)])
        combined[i] = dict(Ne)
    with open(output, "w") as f:
        json.dump(combined, f)
def summarize_SMCpp_results(smcpp_results, output):
    # List to store found files
    combined = {}
    # Walk through the directory and subdirectories
    for i in smcpp_results:
        f=open(i,"r")
        title=f.readline()
        line=f.readline()
        result={}
        while line != "":
            li=line.split(",")
            if li[0] not in result:
                result[li[0]]={li[1]:li[2]}
            else:
                result[li[0]][li[1]] = li[2]
            line=f.readline()

        combined[i]=result

    with open(output, 'w') as json_file:
            json.dump(combined, json_file, indent=4)

def summarize_Stairway_results(stairway_results, generation_interval, output):
    combined={}
    for file in stairway_results:
        Ne = {}  # Store Ne estimates for this file
        with open(file, "r") as f:
            for line in f:
                if line.startswith("mutation_per_site"):  # Skip header
                    continue
                values = line.strip().split("\t")
                try:
                    T = float(values[5])/generation_interval  # 'year' column
                    N = float(values[6])  # 'Ne_median' column
                    N1 = float(values[7])  # Ne_2.5% column
                    N2 = float(values[8])  # Ne_97.5%  column
                    N3 = float(values[9])  # Ne_12.5%	column
                    N4 = float(values[10])  # Ne_87.5%  column

                    Ne.update({T: [N, N1, N2, N3, N4]})
                except ValueError:
                    print(f"Skipping invalid line in {file}: {line}")          
        combined[file]=Ne
    with open(output, 'w') as json_file:
        json.dump(combined, json_file, indent=4)

def summarize_Phlash_results(json_files, output):
    results={}
    for pop_file in json_files:
        with open(pop_file, 'r') as f:
            pop = json.load(f)
        results[pop_file] = pop
    with open(output, 'w') as json_file:
        json.dump(results, json_file, indent=4)

def read_psmc(path, mutation_rate, bin_size):
    """read .psmc, return {T_k(generation): N_k};file is not available using ValueError"""
    with open(path) as f:
        text = f.read()
    if not text:
        raise ValueError("empty or incomplete PSMC output")
    pieces = re.split(r"^//[ \t\r]*$", text, flags=re.M)[:-1]
    blocks = [p for p in pieces if re.search(r"^RS\t", p, flags=re.M)]
    if len(blocks) == 0:
        raise ValueError("no RS lines (empty or incomplete PSMC output)")
    lines = blocks[-1].splitlines()
    tr = [l.split() for l in lines if l.startswith("TR")]
    if len(tr) == 0:
        raise ValueError("no TR line in the last iteration")
    theta = float(tr[0][1])                      # float() recognize 1e-05、nan、inf
    if not math.isfinite(theta) or theta <= 0:
        raise ValueError(f"invalid theta ({tr[0][1]})")
    N_0 = theta / (4 * mutation_rate) / bin_size  # cal per file, avoiding use last file's value

    Ne = {}
    for l in lines:
        if l.startswith("RS"):
            fields = l.split()
            t_k, lambda_k = float(fields[2]), float(fields[3])
            if not (math.isfinite(t_k) and math.isfinite(lambda_k)):
                raise ValueError(f"non-finite value in RS line: {l.strip()}")
            Ne[2 * N_0 * t_k] = N_0 * lambda_k
    return Ne

def summarize_PSMC_results(inputfiles,group_file,mutation_rate,bin_size,output):
    combined={}
    grouped = {}
    popfile_df = pd.read_csv(group_file, sep='\t')
    pop_id = [[f"{row['population']}_{row['group']}", f"{row['sample_id']}.psmc"] for _, row in popfile_df.iterrows()]
    
    for group, element in pop_id:
        if group not in grouped:
            grouped[group] = [element]
        else:
            grouped[group].append(element)

    combined, n_total, n_skipped = {}, 0, 0
    for group, files in grouped.items():
        combined_Ne = {}
        for file in files:
            n_total += 1
            try:
                combined_Ne[file] = read_psmc(file, mutation_rate, bin_size)
            except (ValueError, IndexError, OSError) as e:
                n_skipped += 1
                print(f"WARNING: skipping {file}: {e}", file=sys.stderr)
        if combined_Ne:
            combined[group] = combined_Ne
        else:
            print(f"WARNING: group '{group}' has no valid PSMC files", file=sys.stderr)

    if not combined:
        sys.exit("ERROR: no valid PSMC files in any group")
    print(f"PSMC summary: used {n_total - n_skipped}/{n_total} files", file=sys.stderr)

    with open(output, "w") as json_file:
        json.dump(combined, json_file, indent=4)
def summarize_Relate_results(input_files, output):
    # List to store found files
    combined = {}

    for filename in input_files:
        with open(filename, 'r') as f:
            content = f.readlines()
        groups = re.split(r"\s+", content[0].strip())
        t = re.split(r"\s+", content[1].strip())
        results={}
        for i in range(2,len(content)):
            pop_coal = re.split(r"\s+", content[i].strip())
            group_str = f"{groups[int(pop_coal[0])]} - {groups[int(pop_coal[1])]}"
            if pop_coal[0] == pop_coal[1]:
                c= []
                for j in range(2,len(pop_coal)):
                    if pop_coal[j] == '-nan' or pop_coal[j] == 'nan':
                        c.append(pop_coal[j])
                    else:
                        c.append(0.5/(float(pop_coal[j])+1e-10) )
            else:
                c = pop_coal[2:]
            results[group_str]=dict(zip(t,c))
        
        combined[filename]=results

    with open(output, 'w') as json_file:
        json.dump(combined, json_file, indent=4)

argparser = argparse.ArgumentParser()
argparser.add_argument('-m','--method', required=True)
argparser.add_argument('--mutation', type=float, default=1.25e-8)
argparser.add_argument('--generation_interval', type=float, default=15)
argparser.add_argument('--bin_size', type=float, default=100)
argparser.add_argument('-i','--input_files', nargs='+', required=True)
argparser.add_argument('-g','--group', help="secific to PSMC,group list corresponding to input file")
argparser.add_argument('-o','--output', required=True, help='output file')
args = argparser.parse_args()
input_files = args.input_files
method = args.method
if method == 'GONE':
    summarize_GONE_results(input_files, args.output)
elif method == 'MSMC':
    summarize_MSMC_results(input_files, args.output, args.mutation)
elif method == 'SMCpp':
    summarize_SMCpp_results(input_files, args.output)
elif method == 'stairway':
    summarize_Stairway_results(input_files, args.generation_interval, args.output)
elif method == 'Phlash':
    summarize_Phlash_results(input_files, args.output)
elif method == 'PSMC':
    summarize_PSMC_results(input_files, args.group, args.mutation,args.bin_size,args.output) 
elif method == 'Relate':
    summarize_Relate_results(input_files, args.output)


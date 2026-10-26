#!/usr/bin/env python
"""
Compute Ne(t) from one or more dated tree sequences (.trees.dated.tsz).

For each population in each tree sequence the script computes pairwise
coalescence rates over a shared time grid (either loaded from a file or
built from quantiles of the node-time distribution).  Results are
combined across chromosomes using a weighted average (weights = chromosome
length) and written to a JSON file.

Output JSON structure:
    { "<population_name>": {
        "Weighted_coalrate": { "<t>": <Ne>, ... },
        "<chr_name>":        { "<t>": <Ne>, ... },
        ...
      }, ...
    }

Usage:
    python ts_summary_Ne.py --input *.trees.dated.tsz \
        --timepoint 64 --output tskit_results.json
"""

import argparse
import tskit
import numpy as np
import json
import os
import tszip

parser = argparse.ArgumentParser(description="Summarise Ne from dated tree sequences.")

# Add arguments pathlib.Path
parser.add_argument('--input',     type=str,  nargs="+", required=True, help='Path to the tree sequence')
parser.add_argument('--timepoint', type=str, default="64", help='if no path to the time point, use this value as the number of time points')
parser.add_argument('--output', type=str, required=True, help='Path to the output file')

# Parse arguments
args = parser.parse_args()
# Path to the folder
file_list = args.input
timep=args.timepoint
# Print the list of .trees files
if os.path.exists(timep):
    tq = np.loadtxt(timep, delimiter=",")
    num_time_point = None
elif timep.isdigit():
    num_time_point = int(timep)
    if num_time_point < 1:
        raise ValueError("The number of time points must be a positive integer.")
else:
    print(f"WARNING: Could not recoginize {timep}, defaulting time_point_number=64")
    num_time_point = 64

num_contig = len(file_list)
chr_length = np.zeros(num_contig)
chr_name = []
coalrate_pop = {}
tq_list= []

for count, file in enumerate(file_list):
    chr_name.append(os.path.basename(file).split(".trees.dated.tsz")[0])
    ts = tszip.load(file)
    chr_length[count] = ts.get_sequence_length()
    if num_time_point is not None:
        t = ts.nodes_time
        tq = np.append(0,np.quantile(a=t[t>0], q=[np.linspace(0,1,num_time_point-1)]))
        tq = np.append(tq , np.inf)

    tq_list.append(tq[:-1])

    for pop in ts.populations():
        pop_name = pop.metadata.decode("utf-8")
        pop_name = json.loads(pop_name)['name']
        coal_rate = ts.pair_coalescence_rates(sample_sets=[ts.samples(population=pop.id)],time_windows=tq)
        if pop_name not in coalrate_pop:
            coalrate_pop[pop_name] = [coal_rate]
        else:
            coalrate_pop[pop_name].append(coal_rate)
    
dic_re={}

for key in coalrate_pop:
    coalrate_np = np.array(coalrate_pop[key])
    coalrate_np = np.dot(chr_length, coalrate_np) / np.sum(chr_length)
    tq_np = np.array(tq_list)
    dic_tmp = {}
    dic_tmp["Weighted_coalrate"]=dict(zip(np.mean(tq_np,axis=0),0.5/(coalrate_np+1e-10)))
    for i in range(num_contig):
        dic_tmp[chr_name[i]]=dict(zip(tq_list[i],0.5/(coalrate_pop[key][i]+1e-10)))
    dic_re[key] = dic_tmp

with open(args.output, 'w') as json_file:
    json.dump(dic_re, json_file, indent=4)
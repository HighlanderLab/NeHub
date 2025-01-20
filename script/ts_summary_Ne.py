import argparse
import tskit
import numpy as np
import os
import pathlib
import json
# Create an ArgumentParser object
parser = argparse.ArgumentParser(description="Process some integers.")

# Add arguments pathlib.Path
parser.add_argument('--input',     type=str,  nargs="+", required=True, help='Path to the tree sequence')
parser.add_argument('--timepoint', type=str, default="", help='Path to the time point')
parser.add_argument('--output', type=str, required=True, help='Path to the output file')

# Parse arguments
args = parser.parse_args()

# Path to the folder
file_list = args.input
dic_re={}
timep=args.timepoint
# Print the list of .trees files

for file in file_list:
    ts=tskit.load(file) 
    if timep == "":
        t = ts.nodes_time
        tq = np.append(0,np.quantile(a=t[t>0], q=[np.linspace(0,1,24)]))
        tq = np.append(tq , np.inf)
    else:
        # np.savetxt("timepoint.txt", timepoint, delimiter=",")
        tq = np.loadtxt(timep, delimiter=",")

    coal_rate = ts.pair_coalescence_rates(time_windows=tq)
    dic_re[file]=dict(zip(tq[:-1],0.5/coal_rate))
with open(args.output, 'w') as json_file:
    json.dump(dic_re, json_file, indent=4)


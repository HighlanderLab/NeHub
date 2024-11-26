import argparse
import tskit
import numpy 
# need tskit 0.6.0
import os
import glob

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description="Process some integers.")

# Add arguments
parser.add_argument('--input', type=str, required=True, help='Path to the tree sequence')
parser.add_argument('--timepoint', type=str,default="", help='Path to the time point')
parser.add_argument('--output', type=str, required=True, help='Path to the output file')

# Parse arguments
args = parser.parse_args()

# # Path to the folder
# folder_path = args.input

# # Get all files with the suffix .trees in the folder
# trees_files = glob.glob(os.path.join(folder_path, "*.trees"))

# os.chdir(args.input)

dic_re={}
# Print the list of .trees files
for file in  args.input:

    ts=tskit.load(file) 
    if args.timepoint == "":
        t = ts.nodes_time
        tq = np.append(0,np.quantile(a=t[t>0], q=[np.linspace(0,1,24)]))
        tq = np.append(tq , np.inf)
        coal_rate = ts.pair_coalescence_rates(time_windows=tq)
    else:
        # np.savetxt("timepoint.txt", timepoint, delimiter=",")
        tq = np.loadtxt(args.timepoint, delimiter=",")
        coal_rate = ts.pair_coalescence_rates(time_windows=tq)
    dic_re[file]=dict(zip(tq[:-1],0.5/coal_rate))

with open(arg.output, 'w') as json_file:
    json.dump(dic_re, json_file, indent=4)
import os
import pickle
import argparse
import numpy as np
import json
import glob

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description="Process some integers.")

# Add arguments
parser.add_argument('--input', type=str,nargs="+",  required=True, help='Path to the result folder')
parser.add_argument('--output',type=str, required=True, help='output folder')
# Parse arguments
args = parser.parse_args()
combined_li = {}
for sampling_dir in args.input:
    try:
        # Get all .pkl files in the directory
        pkl_files = glob.glob(os.path.join(sampling_dir, '*.pkl'))
    except FileNotFoundError:
        print(f"Directory not found: {sampling_dir}")
        break
    results={}
    for pop_file in pkl_files:
        try:
            # Load the pickle file
            with open(pop_file, 'rb') as f:
                pop_results = pickle.load(f)

            times = np.array([dm.eta.t[1:] for dm in pop_results])
            # choose a grid of points at which to evaluate the size history functions
            T = np.geomspace(times.min(), times.max(), 1000)
            Nes = np.array([dm.eta(T, Ne=True) for dm in pop_results])
            results[os.path.basename(pop_file)]=dict(zip(T,np.median(Nes, axis=0)))
        except Exception as e:
            print(f"Error processing {pop_file}: {e}")
            continue

    combined_li[sampling_dir] = results

with open(args.output, 'w') as json_file:
    json.dump(combined_li, json_file, indent=4)


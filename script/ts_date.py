import tsdate
import tskit
import argparse
import os
import tsinfer
# Create an ArgumentParser object
parser = argparse.ArgumentParser(description="Process some integers.")

# Add arguments
parser.add_argument('--input', type=str, required=True, help='Path to the sampledata file')
parser.add_argument('--mutation_rate', type=float, required=True, help='mutation rate')
# parser.add_argument('--output', type=str, required=True, help='Path to the output trees file')
# Parse arguments
args = parser.parse_args()

input_tree=args.input
mutation = args.mutation_rate
output = input_tree.replace(".trees",".dated.trees")

ts=tskit.load(input_tree)
ts=ts.simplify()
dated_ts = tsdate.date(ts, mutation_rate=mutation)
dated_ts.dump(output)
import tsdate
import tskit
import argparse
import os
# Create an ArgumentParser object
parser = argparse.ArgumentParser(description="Process some integers.")

# Add arguments
parser.add_argument('--input', type=str, required=True, help='Path to the sampledata file')
parser.add_argument('--mutation_rate', type=str, required=True, help='mutation rate')
# parser.add_argument('--output', type=str, required=True, help='Path to the output trees file')
# Parse arguments
args = parser.parse_args()

samples=args.input
mutation = args.mutation_rate
output = samples.replace(".samples",".dated.trees",)
ts = tsinfer.infer(samples)

print("Inferred tree sequence: {} trees over {} Mb ({} edges)".format(
        ts.num_trees, ts.sequence_length / 1e6, ts.num_edges
                ))

redated_ts = tsdate.date(ts, mutation_rate=mutation)
redated_ts.dump(output)
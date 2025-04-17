import tsdate
import tsinfer

import argparse

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description="Process some integers.")

# Add arguments
parser.add_argument('--samples', type=str, required=True, help='Path to the sample file from generate_sample_file.py')
parser.add_argument('--output', type=str, required=True, help='Path to the output file (.trees)')
parser.add_argument("--mutation", type=float, default=1e-08, help='mutation rate')
# Parse arguments
args = parser.parse_args()


# Do the inference
ts = tsinfer.infer(args.samples)

print(
    "Inferred tree sequence: {} trees over {} Mb ({} edges)".format(
        ts.num_trees, ts.sequence_length / 1e6, ts.num_edges
    )
)

redated_ts = tsdate.date(sim_ts, mutation_rate=args.mutation)

ts.dump(args.output)
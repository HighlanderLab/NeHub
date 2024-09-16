import random
import argparse
import os


parser = argparse.ArgumentParser(description="Process input data and write results to output file")

# Add arguments for input and output file paths
parser.add_argument("-i", "--input", required=True, help="Path to the input data file")
parser.add_argument("-o", "--output", required=True, help="Path to the output file where results will be written")
parser.add_argument("-n", "--num_sampling", type=int, default=10, help="Number of bootstrap samples to generate")
parser.add_argument("-s", "--sample_size", type=int, default=10, help="Number of elements in each bootstrap sample")
# Parse the command-line arguments
args = parser.parse_args()
# Define the total set of numbers (for example, from 1 to 40)
li=[]
with open(args.input, "r") as f:
    line = f.readline()
    while line!="":
        li.appand(line)
        line = f.readline()


def bootstrap_sample(data, sample_size, num_sampling):
    # Generate 'num_samples' bootstrap samples
    bootstrap_samples = []
    for i in range(num_sampling):
        sample = [random.choice(data) for _ in range(sample_size)]
        bootstrap_samples.append(sample)
    return bootstrap_samples


# Perform bootstrap sampling
num_sampling = args.num_sampling  # Number of bootstrap samples
sample_size = args.sample_size  # Number of elements in each sample
bootstrap_samples = bootstrap_sample(li, sample_size, num_sampling)

# Write the results to the output file
for i, sample in enumerate(bootstrap_samples):
    with open(os.path.join(args.output,f"sampling{i}.txt"), "w") as f:
        
        f.write(sample.join(""))
import random
import argparse
import os


parser = argparse.ArgumentParser(description="Process input data and write results to output file")

# Add arguments for input and output file paths
parser.add_argument("-i", "--input", required=True, help="Path to the input data file")
parser.add_argument("-o", "--output", required=True, help="Path to the output file where results will be written")
parser.add_argument("-n", "--num_sampling", type=int, default=10, help="Number of bootstrap samples to generate")
parser.add_argument("-s", "--sample_size", type=int, default=10, help="Number of elements in each bootstrap sample")
parser.add_argument("-m", "--sampling_method", type=str, default="bootstrap", help="There is two method to resampling, bootstraps and jackknife")
parser.add_argument("-b", "--block_size", type=int, default=1, help="Block size, default value is 1")


# Parse the command-line arguments
args = parser.parse_args()
# Define the total set of numbers (for example, from 1 to 40)
dic={}
with open(args.input, "r") as f:
    title = f.readline()
    line = f.readline()
    while line!="":
        ind=line.strip("\n").split("\t")[0]
        pop=line.strip("\n").split("\t")[1]
        if ind in dic:
            print("Please check id, there is repeat individual")
            break
        else:
            dic[ind]=pop
        line = f.readline()
        


def bootstrap_sample(data, sample_size, num_sampling):
    # "boostrap without replacement"
    # Generate 'num_samples' bootstrap samples
    bootstrap_samples = []
    for i in range(num_sampling):
        sample = random.sample(data, sample_size)  
        bootstrap_samples.append(sample)
    return bootstrap_samples

def Jackknife_sample(data, block_size):
    jackknife_samples = []
    n = len(data)
    if block_size <= 0 or block_size > n:
        raise ValueError("Block size must be between 1 and the length of the dataset.")
    block_jackknife_samples = []
    num_blocks = n // block_size
    # Generate a jackknife sample by leaving out one element at a time
    for i in range(num_blocks):
        start = i * block_size
        end = start + block_size
        sample = data[:start] + data[end:]  # Leave out the block from start to end
        block_jackknife_samples.append(sample)
    
    return block_jackknife_samples

if args.sampling_method == "bootstrap":

    # Perform bootstrap sampling
    num_sampling = args.num_sampling  # Number of bootstrap samples
    sample_size = args.sample_size  # Number of elements in each sample
    li=[key for key in dic.keys()]
    bootstrap_samples = bootstrap_sample(li, sample_size, num_sampling)

    # Write the results to the output file
    for i, sample in enumerate(bootstrap_samples):
        with open(os.path.join(args.output,f"sampling{i}.popLabel"), "w") as f:
            f.write(title)
            for k in sample:
                f.write(f"{k}\t{dic[k]}\n")
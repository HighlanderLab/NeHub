import random
import argparse
import os
import cyvcf2

def readvcf(vcf,sample):
    dic_output={}
    pos=0
    for record in vcf:  # Loop over variants, each assumed at a unique site
        if pos==0:
            i_index = [vcf.samples.index(s) for s in sample]
        
        chr=record.CHROM
        position = record.POS
        interval = position-pos
        pos=position

        variants_alleles=""
        variants_alleles2=""
        for alleles in record.gt_bases[i_index]:
            variants_alleles = variants_alleles+alleles[0]+alleles[2]
            if "|" not in alleles:
                variants_alleles2 = variants_alleles2+alleles[2]+alleles[1]                 
            else:
                variants_alleles2 = variants_alleles

        if variants_alleles ==variants_alleles2:
            final_col = variants_alleles
        else:
            final_col = variants_alleles +"," + variants_alleles2
        if chr not in dic_output:
            dic_output[chr]=[f"{chr}\t{position}\t{interval}\t{final_col}\n"]
        else:
            dic_output[chr].append(f"{chr}\t{position}\t{interval}\t{final_col}\n")
    return dic_output

parser = argparse.ArgumentParser(description="Process input data and write results to output file")

# Add arguments for input and output file paths
parser.add_argument("-i", "--input", required=True, help="Path to the input data file, vcf file")
parser.add_argument("-p", "--popLabel", required=True, help="Path to the input data file, pop label file")
parser.add_argument("-w", "--workDir", required=True, help="Path to the output file where results will be written, tab for msmc")
parser.add_argument("-n", "--num_sampling", type=int, default=5, help="Number of bootstrap samples to generate")
parser.add_argument("-s", "--sample_size", type=int, default=4, help="Number of elements in each bootstrap sample,recommondate 2,4,8")


# Parse the command-line arguments
args = parser.parse_args()
# Define the total set of numbers (for example, from 1 to 40)
dic={}
with open(args.popLabel, "r") as f:
    title = f.readline()
    line = f.readline()
    while line!="":
        ind=line.strip("\n").split("\t")[0]
        pop=line.strip("\n").split("\t")[1]
        if pop in dic:
            dic[pop].append(ind)
        else:
            dic[pop]=[ind]
        line = f.readline()

for i in range(args.num_sampling):   
    # Check if the folder exists
    if not os.path.exists(os.path.join(args.workDir,f"./sampling{i}/")):
        # Create the folder
        os.mkdir(os.path.join(args.workDir,f"./sampling{i}/"))
    f = open(os.path.join(args.workDir,f"./sampling{i}/pop.Label"), "w") 
    f.write(title)
    vcf = cyvcf2.VCF(args.input)

    count=0
    parameters = []
    for pop in dic:
        sample = random.sample(dic[pop], args.sample_size) 
        for k in sample:
            f.write(f"{k}\t{pop}\n")
        f.close()
        dic_chr=readvcf(vcf,sample)
        for chr in dic_chr:
            with open(os.path.join(args.workDir,f"./sampling{i}/msmc_{chr}.tab"), "w") as f2:
                f2.write("".join(dic_chr[chr]))
        
        parameters.append([str(i+count*args.sample_size*2) for i in range(args.sample_size*2)])
        count+=1

with open(os.path.join(args.workDir,"finish_flag"), "w") as f3:
    for i in parameters:
        f3.write(",".join(i))
        f3.write("\n")
    if count>1:
        cross_coal=[]
        for j in parameters[0]:
            for k in parameters[1]:
                cross_coal.append(j+"-"+k)
        f3.write(",".join(cross_coal))


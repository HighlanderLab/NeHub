# transfer file to vcf format
import tskit
import argparse

parser = argparse.ArgumentParser(description="Process input data and write results to output file")

# Add arguments for input and output file paths
parser.add_argument("-i", "--input", required=True, help="Path to the input data file")
parser.add_argument("-r", "--recomb",type=float,default=9.5e-9,help="recombination rate")
# Parse the command-line arguments
args = parser.parse_args()

ts=tskit.load(args.input)
length=ts.sequence_length
n_dip_indv = int(ts.num_samples / 2)
indv_names = [f"tsk_{i}indv" for i in range(n_dip_indv)]

import gzip

with gzip.open(f"{args.input}.vcf.gz", "wt") as vcf_file:
    ts.write_vcf(vcf_file, individual_names=indv_names,contig_id="chr1")

import random
def reference_genome(ts,length,output):
    s=[random.choices(["A","T","C","G"],[0.225,0.225,0.275,0.275],k=1)[0] for i in range(length)]
    for i in ts.sites():
        s[int(i.position)-1]=i.ancestral_state
    # write reference genome, non-variant sites will be written randomly (ATCG)
    with open(output,"w") as f:
        f.write(">chr1\n")
        for i in range(len(s)):
            if (i%60==0) and (i!=0):
                f.write("\n")
            f.write(s[i])
reference_genome(ts,int(length),f"{args.input}.fa")

with open(f"{args.input}.popLabel","w") as f:
    f.write("sample_id\tpopulation\n")
    for i in ts.individuals():
        f.write("{}\t{}\n".format(f"tsk_{i.id}indv",i.population))

r=args.recomb*1e+8
lp=0
gm=0
with open(f"{args.input}.genetic_map.txt","w") as f:
    f.write("pos\tCOMBINED_rate\tGenetic_Map\n")
    
    for i in ts.sites():
        p=int(i.position)
        gm += (p-lp)*r/1e+6
        lp=int(i.position)
        f.write(f"{lp}\t{r}\t{gm}\n")
        
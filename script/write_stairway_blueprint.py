blueprint_config=[
"#example blueprint file",
"#input setting",
"popid: two-epoch # id of the population (no white space)",
"nseq: 30 # number of sequences",
"L: 10000000 # total number of observed nucleic sites, including polymorphic and monomorphic",
"whether_folded: false # whethr the SFS is folded (true or false)",
"SFS: 9513.26	3796.47	2106.24	1351.505	962.91	736.865	597.75	499.845	429.99	381.62	341.965	310.27	284.87	261.28	242.11	230.055	217.975	204.48	197.28	185.545	176.055	167.62	161.55	154.935	147.83	142.17	137.31	133.3	124.955 # snp frequency spectrum: number of singleton, number of doubleton, etc. (separated by white space)",
"#smallest_size_of_SFS_bin_used_for_estimation: 1 # default is 1; to ignore singletons, uncomment this line and change this number to 2",
"#largest_size_of_SFS_bin_used_for_estimation: 29 # default is n-1; to ignore singletons, uncomment this line and change this number to nseq-2",
"pct_training: 0.67 # percentage of sites for training",
"nrand: 7	15	22	28 # number of random break points for each try (separated by white space)",
"project_dir: two-epoch # project directory",
"stairway_plot_dir: stairway_plot_es # directory to the stairway plot files",
"ninput: 200 # number of input files to be created for each estimation",
"#random_seed: 6",
"#output setting",
"mu: 1.2e-8 # assumed mutation rate per site per generation",
"year_per_generation: 24 # assumed generation time (in years)",
"#plot setting",
"plot_title: two-epoch # title of the plot",
'xrange: 0.1,10000 # Time (1k year) range; format: xmin,xmax; "0,0" for default',
'yrange: 0,0 # Ne (1k individual) range; format: xmin,xmax; "0,0" for default',
"xspacing: 2 # X axis spacing",
"yspacing: 2 # Y axis spacing",
"fontsize: 12 # Font size"]
import argparse
import re
import os 
# Create argument parser
parser = argparse.ArgumentParser(description="build blueprint for stairway plot")

# Add arguments
parser.add_argument("-i", "--input", type=str, required=True, help="Input file path")
parser.add_argument("-stairway", "--stairway_plot_dir", type=str, required=True, help="The path of stairway_plot_dir")
parser.add_argument("-o", "--output", type=str, required=True, help="Output file path")
parser.add_argument("--nind", type=int, required=True, help="the number of individual")
parser.add_argument("--length", type=int,required=True, help="total number of observed nucleic sites, including polymorphic and monomorphic")
parser.add_argument("--chr", type=str,required=True, help="chrmosome name")
parser.add_argument("--mutation_rate", type=float,required=True, help="mutation rate")
parser.add_argument("--fold", type=int,default=1, help="Whether SFS is generated from fold. Default 1 is unfold.")
parser.add_argument("--ploid", type=int,default=2, help="Ploidy level (e.g. 2 for diploid humans)")
parser.add_argument("--generation_time", type=float,default=24, help="generation time, # summary results use 'generation' as units, so this value didn't affect final results, but if you are insterested the plot, it must be setting")
# Parse arguments
args = parser.parse_args()

with open(args.input,'r') as f:
    lines = f.readlines()
li = lines[2].strip("\n").split(" ")
nseq=args.ploid*args.nind
fold="true"
if args.fold == 1:
    sfs="\t".join(li[1:args.nind+2])
else:
    sfs="\t".join(li[1:-1])
    fold="false"

inputname=os.path.basename(args.input)
inputname = inputname.split("_SFS.obs")[0]

# blueprint_config[2]= re.sub(r":(.*?)#", f": {inputname[1]} #", blueprint_config[2])
# because of limits of snakemake, now use the same name 
blueprint_config[2]= re.sub(r":(.*?)#", f": results #", blueprint_config[2])
blueprint_config[3]= re.sub(r":(.*?)#", f": {nseq} #", blueprint_config[3])
blueprint_config[4]= re.sub(r":(.*?)#", f": {args.length} #", blueprint_config[4])
blueprint_config[5]= re.sub(r":(.*?)#", f": {fold} #", blueprint_config[5])
blueprint_config[6]= re.sub(r":(.*?)#", f": {sfs} #", blueprint_config[6])
blueprint_config[11]= re.sub(r":(.*?)#", f": {inputname} #", blueprint_config[11])
blueprint_config[12]= re.sub(r":(.*?)#", f": {args.stairway_plot_dir} #", blueprint_config[12])
blueprint_config[16]= re.sub(r":(.*?)#", f": {args.mutation_rate} #", blueprint_config[16])
blueprint_config[17]= re.sub(r":(.*?)#", f": {args.generation_time} #", blueprint_config[17])
blueprint_config[19]= re.sub(r":(.*?)#", f": {inputname} #", blueprint_config[19])
print(inputname)
with open(args.output,"w") as f:
    f.write("\n".join(blueprint_config))
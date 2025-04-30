
input_parameters_file=[
'#INPUT_PARAMETERS_FILE is only used to GONE',
'',
'########################################################',
'',
'PHASE=2 ### Phase = 0 (pseudohaploids), 1 (known phase), 2 (unknown phase)',
'cMMb=1  ### CentiMorgans per Megabase (if distance is not available in map file).',
'DIST=1  ### none (0), Haldane correction (1) or Kosambi correction (2)',
'NGEN=2000 ### Number of generations for which linkage data is obtained in bins',
'NBIN=400  ### Number of bins (e.g. if 400, each bin includes NGEN/NBIN = 2000/400 = 5 generations)',
'MAF=0.0   ### Minor allele frequency (0-1) (recommended 0)',
'ZERO=1    ### 0: Remove SNPs with zeroes (1: allow for them)',
'maxNCHROM=-99  ### Maximum number of chromosomes to be analysed (-99 = all chromosomes; maximum number is 200)',
'maxNSNP=50000 ### Maximum approx number of SNPs per chromosomes to be analysed (maximum number is 50000)',
'hc=0.05   ### Maximum value of c analysed (recommended 0.05; maximum is 0.5)',
'REPS=40   ### Number of replicates to run GONE (recommended 40)',
'threads=-99  ### Number of threads (if -99 it uses all possible processors)',
'',
'###################################################################'
]

import argparse
import re
import os 
# Create argument parser
parser = argparse.ArgumentParser(description="build input parameters file for GONE")

# Add arguments
parser.add_argument("-o", "--output", type=str, required=True, help="Output file path")
parser.add_argument("-p", "--phase", type=str, required=True, help="whether the phase is known or not, 0: pseudohaploids, 1: known phase, 2: unknown phase")
# Parse arguments
args = parser.parse_args()
input_parameters_file[4]= re.sub(r"=(.*?)###", f"={args.phase} ###", input_parameters_file[4])
# we don't change any value, because this setting set is suitable for most of analysis. 
# if you are interested into change value we give the template to modify. 

'''
# before parser.parse_args() insert :
parser.add_argument("--chr", type=str,help="chrmosome name")
# after parser.parse_args() insert:
input_parameters_file[2]= re.sub(r"=(.*?)###", f": {args.xxx} #", input_parameters_file[2])
'''

with open(args.output,"w") as f:
    f.write("\n".join(input_parameters_file))



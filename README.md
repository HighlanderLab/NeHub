# ne-traj-pipe
Pipeline to infer population size trajectories using a variety of tools

Sampling


|Tools| input| link|
|-----|-----|---|
|GONE | plink file .ped .bim|https://github.com/esrud/GONE|
|PSMC|fastq |https://github.com/lh3/psmc|
|SMC++ | vcf.gz|https://github.com/popgenmethods/smcpp|
|pHlash | vcf.gz| https://github.com/jthlab/phlash|
|Relate (ARG)| ts or (anc.gz and mut.gz)|https://myersgroup.github.io/relate/|
|IBDNe |plink file .bim .bed .fam|https://github.com/hennlab/AS-IBDNe|
|Fastsimcoal2|sfs|https://cmpg.unibe.ch/software/fastsimcoal2/|
|Dadi |sfs|https://dadi.readthedocs.io/en/latest/|

# Run on Eddie, interactive node
```
module load anaconda
conda activate HLab_tsinfer
qlogin -l h_vmem=4G
cd ./ne-traj-pipe
snakemake -j 4 -s snakefile
```

# If you need reset one of input in config
```
snakemake -j 4 -s snakefile --config vcffile="dec0.vcf.gz"
# config in command line is proior, so it will cover the original one.
```

# cluster
before run cluster, please check drmaa path in your ```~/.bash_profile``` or equivalent: ```export DRMAA_LIBRARY_PATH=/exports/applications/gridengine/ge-8.6.5/lib/lx-amd64/libdrmaa.so```
```
snakemake -j 8 -profile profile --config vcffile="dec0.vcf.gz" &
```
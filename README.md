# ne-traj-pipe
Pipeline to infer population size trajectories using a variety of tools

Sampling


|Tools| input| link| install approach|
|-----|-----|---|---|
|GONE | plink file .ped .bim|https://github.com/esrud/GONE|git clone|
|PSMC|fastq |https://github.com/lh3/psmc|conda install|
|MSMC2|fastq |https://github.com/stschiff/msmc2|conda install|
|SMC++ | vcf.gz|https://github.com/popgenmethods/smcpp|conda install|
|pHlash | vcf.gz| https://github.com/jthlab/phlash|pip install|
|tskit(tsinfer and tsdate)|ts|https://tskit.dev/|conda install|
|Relate (ARG)| haps and sample|https://myersgroup.github.io/relate/|package|
|IBDNe |plink file .bim .bed .fam|https://github.com/hennlab/AS-IBDNe|java package|
|Fastsimcoal2|sfs|https://cmpg.unibe.ch/software/fastsimcoal2/|conda install|
|Dadi |sfs|https://dadi.readthedocs.io/en/latest/|conda install|
|stairway-plot-v2|sfs|https://github.com/xiaoming-liu/stairway-plot-v2/tree/master|package|

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
Also must login a node, otherwise report error 
```
OpenBLAS blas_thread_init: pthread_create failed for thread 7 of 16: Resource temporarily unavailable
OpenBLAS blas_thread_init: RLIMIT_NPROC 200 current, 200 max
```

```
snakemake --profile profile --config vcffile="./test_data/dec0" &
```
OR to run snakemake in the background, with output redirected to a file, and prevent the process from hanging after the command is sent to the background
```
snakemake --profile profile --config vcffile="./test_data/dec0" > ../script/GONE_pipeline/snakemake_output.log 2>&1 &
```

## check the currently running process
```
ps u -u [username]
```

## Inference environments on Eddie
### PSMC
`source activate ne_pipe_psmc`  
Origianally created as:  
`conda create -n ne_pipe_psmc -c bioconda -c conda-forge psmc`

### MSMC2
`source activate ne_pipe_msmc2`  
Originally created as:  
`conda create -n ne_pipe_msmc2 -c bioconda -c conda-forge msmc2`

### SMC++
Should be run via singularity
```
module load singularity/3.11.4 # or similar
export SINGULARITY_TMPDIR="\$TMPDIR"
singularity run docker://terhorst/smcpp:latest -h # replace -h by SMC++ parameters
```

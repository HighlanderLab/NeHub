# ne-traj-pipe
Pipeline to infer population size trajectories using a variety of tools

Sampling samples

|Tools| input| link| install approach|if built-in boostrap for chunk of chromosome|
|-----|-----|---|---|----|
|GONE | plink file .ped .bim|https://github.com/esrud/GONE|git clone|yes, but sampling variants independly|
|PSMC|fastq |https://github.com/lh3/psmc|conda install|yes|
|MSMC2|vcf |https://github.com/stschiff/msmc2; https://github.com/stschiff/msmc-tools/ |conda install, git clone|yes|
|SMC++ | vcf.gz|https://github.com/popgenmethods/smcpp|conda install||
|pHlash | vcf.gz| https://github.com/jthlab/phlash|pip install||
|tskit(tsinfer and tsdate)|ts|https://tskit.dev/|conda install||
|Relate (ARG)| haps and sample|https://myersgroup.github.io/relate/|package||
|SINGER |vcf |https://github.com/popgenmethods/SINGER|package ||
|HapNe-LD/HapNe-IBD |vcf|https://github.com/palamaraLab/hapne|pip install||
|Dadi |sfs|https://dadi.readthedocs.io/en/latest/|conda install||
|Fastsimcoal2|sfs|https://cmpg.unibe.ch/software/fastsimcoal2/|conda install||
|InferNo|.dat?|https://github.com/ZhendongHuang/Inference for demographic history|git clone | |

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
## prepare
before run cluster, please check drmaa path in your ```~/.bash_profile``` or equivalent: ```export DRMAA_LIBRARY_PATH=/exports/applications/gridengine/ge-8.6.5/lib/lx-amd64/libdrmaa.so```
Also must login a node, otherwise report error 
```
OpenBLAS blas_thread_init: pthread_create failed for thread 7 of 16: Resource temporarily unavailable
OpenBLAS blas_thread_init: RLIMIT_NPROC 200 current, 200 max
```

check if conda.sh is included in ```~/.bashrc```, if not it will report error
```/bin/bash: conda: command not found```

## run snakemake

```
snakemake --profile absolute/path/to/profile --config vcffile="./test_data/dec0" &
```
OR to run snakemake in the background, with output redirected to a file, and prevent the process from hanging after the command is sent to the background
```
snakemake --profile profile --config vcffile="./test_data/dec0" > ../script/GONE_pipeline/snakemake_output.log 2>&1 &
```

## check the currently running process
```
qstat
```

## Inference environments on Eddie
### PSMC
`source activate ne_pipe_psmc`  
Origianally created as:  
`conda create -n ne_pipe_psmc -c bioconda -c conda-forge psmc`

### MSMC2 
`source activate ne_pipe_msmc2`  
`msmc2_Linux`
Originally created as:  
`conda create -n ne_pipe_msmc2 -c bioconda -c conda-forge msmc2`


### SMC++ 
Should be run via singularity
```
module load singularity/3.11.4 # or similar
export SINGULARITY_TMPDIR="/exports/eddie/scratch/username/"
singularity run docker://terhorst/smcpp:latest -h # replace -h by SMC++ parameters
```
run:
```
module load singularity
export SINGULARITY_TMPDIR="/home/s2587593/PigData/NeTest"
singularity run docker://terhorst/smcpp:latest vcf2smc ./test_data/con0.vcf.gz ../NeTest/con0/chr1.smc.gz chr1 
```

Error:
```
FATAL:   Unable to handle docker://terhorst/smcpp:latest uri: while building SIF from layers: while creating squashfs: create command failed: exit status 1: 
FATAL ERROR: Out of memory (cache_alloc)
```
Solution:
```
qlogin use 16G memory.
```

### dadi-cli
`conda create -n ne_pipe_dadiCli -c conda-forge dadi-cli` 

Error:
```
/exports/cmvm/eddie/eb/groups/HighlanderLab/anaconda/envs/ne_pipe_dadiCli/bin/python3.12: can't open file '/exports/cmvm/eddie/eb/groups/HighlanderLab/anaconda/envs/ne_pipe_dadiCli/bin/dadi-cli': [Errno 13] Permission denied
```

### fastsimcoal2
`conda create -n ne_pipe_fastsimcoal2 -c bioconda fastsimcoal2`

### phlash
phlash requires a GPU to run fast. THere are GPUs available on Eddie, but we probably have to play around a little to figure out how to use these. Hannes has set up two versions.

#### No GPU
This is not a conda environment! Python venv in:  
`/exports/cmvm/eddie/eb/groups/HighlanderLab/share/venvs/ne_pipe_phlash_nogpu`

To activate, run:  
```
source /exports/cmvm/eddie/eb/groups/HighlanderLab/share/venvs/ne_pipe_phlash_nogpu/bin/activate
# then, try
python
>>> import phlash
```


```
##Setup code (so we remember)
#python -m venv --copies ne_pipe_phlash_nogpu
#source ne_pipe_phlash_nogpu/bin/activate
#python -m pip install -U pip setuptools
#python -m pip install 'phlash'
```

#### With GPU support
This is not a conda environment! Python venv in:  
`/exports/cmvm/eddie/eb/groups/HighlanderLab/share/venvs/ne_pipe_phlash_gpu`
 
To activate, run:  
```
source /exports/cmvm/eddie/eb/groups/HighlanderLab/share/venvs/ne_pipe_phlash_gpu/bin/activate
# then, try
python
>>> import phlash
```

```
##Setup code
#python -m venv --copies ne_pipe_phlash_gpu
#source ne_pipe_phlash_gpu/bin/activate
#python -m pip install -U pip setuptools
#python -m pip install 'phlash[gpu]'
```
#### After source I have this erorr report

```
$ python
Could not find platform independent libraries <prefix>
Could not find platform dependent libraries <exec_prefix>
Python path configuration:
  PYTHONHOME = (not set)
  PYTHONPATH = (not set)
  program name = 'python'
  isolated = 0
  environment = 1
  user site = 1
  safe_path = 0
  import site = 1
  is in build tree = 0
  stdlib dir = '/home/hbecher/miniconda3/lib/python3.11'
  sys._base_executable = '/home/hbecher/miniconda3/bin/python'
  sys.base_prefix = '/home/hbecher/miniconda3'
  sys.base_exec_prefix = '/home/hbecher/miniconda3'
  sys.platlibdir = 'lib'
  sys.executable = '/exports/cmvm/eddie/eb/groups/HighlanderLab/share/venvs/ne_pipe_phlash_gpu/bin/python'
  sys.prefix = '/home/hbecher/miniconda3'
  sys.exec_prefix = '/home/hbecher/miniconda3'
  sys.path = [
    '/home/hbecher/miniconda3/lib/python311.zip',
    '/home/hbecher/miniconda3/lib/python3.11',
    '/home/hbecher/miniconda3/lib/python3.11/lib-dynload',
  ]
Fatal Python error: init_fs_encoding: failed to get the Python codec of the filesystem encoding
Python runtime state: core initialized
ModuleNotFoundError: No module named 'encodings'

Current thread 0x00007ffa86743740 (most recent call first):
  <no Python frame>
```

# For real data, mask used to exclude Non-neutral sites and errors
tutorial from paper (Exploring population size changes using SNP frequency spectra)

1. at least 50 kb away from any coding region
2. mapping errors
3. ancestral allele with low confidence 

4. CpG Islands: CpG-rich regions often exhibit higher mutation rates due to spontaneous deamination. 
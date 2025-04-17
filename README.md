# ne-traj-pipe
A comprehensive pipeline for inferring population size trajectories using multiple demographic inference tools.
## Table of Contents
- [Available Tools](#available-tools)
- [Configuration](#configuration)
- [Data Processing](#data-processing)
- [Usage](#usage)
- [Cluster Usage](#cluster-usage)
- [Installation](#installation)


## Available Tools
|Tools| input| link| install approach|if built-in boostrap for chunk of chromosome|Avaliable now|
|-----|-----|---|---|----|----|
|GONE | plink file .ped .bim|https://github.com/esrud/GONE|git clone|yes, but sampling variants independly|Yes|
|PSMC|fastq |https://github.com/lh3/psmc|conda install|yes|Yes|
|MSMC2|vcf |https://github.com/stschiff/msmc2; https://github.com/stschiff/msmc-tools/ |conda install, git clone|yes|Yes|
|SMC++ | vcf.gz|https://github.com/popgenmethods/smcpp|conda install|No|Yes|
|pHlash | vcf.gz| https://github.com/jthlab/phlash|pip install|No|Yes|
|tskit(tsinfer and tsdate)|ts|https://tskit.dev/|conda install|No|Yes|
|Relate (ARG)| haps and sample|https://myersgroup.github.io/relate/|package|No|Yes|
|Stairway Plot v2|SFS|https://github.com/xiaoming-liu/stairway-plot-v2|No|Yes|
|HapNe-LD/HapNe-IBD |vcf|https://github.com/palamaraLab/hapne|pip install|No|Now it only support human|
|SINGER |vcf |https://github.com/popgenmethods/SINGER|package ||It is for sampling and resampling ARGs|
|Dadi |sfs|https://dadi.readthedocs.io/en/latest/|conda install||-|
|Fastsimcoal2|sfs|https://cmpg.unibe.ch/software/fastsimcoal2/|conda install||-|
|InferNo|.dat?|https://github.com/ZhendongHuang/Inference for demographic history|git clone| |-|

### Memory Requirements
- SMC++ requires at least 16GB memory
- pHlash benefits from GPU acceleration

## Configuration

<a name="global-parameters"></a>
## Global Parameters ([`global_config.yaml`](config/global_config.yaml))

- `vcffile`: Path to the input prefix of VCF file (must end with .vcf.gz, for example "test". don't suggest 
  including "." in file name, prefix file name of popLabel(.popLabel) must be same with vcffile)
- `ref_genome`: Path to the reference genome FASTA file
- `ancestral_genome`: Path to the ancestral genome FASTA file
- `workdir`: Working directory for output files
- `numSampling`: Number of sampling iterations (set to 1 if no sampling needed)
- `sampleSize`: Number of samples to use in each iteration
- `Method`: List of methods to run (e.g., ['GONE', 'PSMC', 'Tskit', 'MSMC', 'SMCpp', 'Relate', 'Phlash', 'Stairway'])
- `mutation_rate`: Mutation rate per site per generation (default: 4.8e-9)

### Global Parameters
The pipeline uses YAML configuration files located in the `config/` directory. The main configuration file is [`global_config.yaml`](config/global_config.yaml):

```yaml
vcffile: "/path/to/your/vcf/prefix"  # Must match popLabel prefix
ref_genome: "/path/to/reference.fa"
ancestral_genome: "/path/to/ancestral.fa"
workdir: "/path/to/output/directory"
numSampling: 5
sampleSize: 80
Method: ['GONE', 'PSMC', 'Stairway']  # Available methods
mutation_rate: 4.8e-9
```

### Tool-Specific Configurations

### Stairway Plot (`Stairway_config.yaml`)
- `easysfs_path`: Path to the easySFS.py executable
  - Used for converting VCF to SFS format
  - Must be the full path to the easySFS script
- `stairway_config`: Path to the Stairway Plot JAR file
  - Contains the Stairway Plot executable and resources
- `fold`: SFS folding option
  - 1: Use folded SFS (default)
  - 0: Use unfolded SFS
  - Folding is recommended when ancestral states are uncertain
- `ploids`: Ploidy level of the organism
  - 2: Diploid (default for most organisms)
  - Other values for polyploid species
- `generation_interval`: Generation time in years
  - Used for converting results to calendar years
  - Does not affect the actual analysis (because the final results are along generation), only the interpretation of results
  - Default: 5 years


#### MSMC (`MSMC_config.yaml`)
- `num_sampling_msmc`: 20
- `sample_size_msmc`: 4
- `timeSegmentPattern`: "1*2+25*1+1*2+1*3"
- `rhoOverMu`: 0.22
- `maxIterative`: 4

#### PSMC (`PSMC_config.yaml`)
- `psmc_setting`: "4+50*1+4+6"
- `bin_size`: 100

#### GONE (`GONE_config.yaml`)
- `GONE_path`: Path to GONE installation
- `GONE_setting`: Path to parameters file

#### Relate (`Relate_config.yaml`)
- `Relate_path`: Path to Relate installation
- `N`: 20000

#### SMC++ (`SMCpp_config.yaml`)
- `contig`: chr1

#### tskit (`tskit_config.yaml`)
- `timepoint`: Time point file

## Data Processing
For real data, mask used to exclude Non-neutral sites and errors
tutorial from paper (Exploring population size changes using SNP frequency spectra)
### Site Filtering Recommendations
For real data analysis, consider masking:
1. Regions within 50kb of coding regions
2. Mapping errors
3. Low-confidence ancestral alleles
4. CpG Islands (due to higher mutation rates, CpG-rich regions often exhibit higher mutation rates due to spontaneous deamination.)

## Usage

If you define all parameters in config, you can run directly.
```bash
module load anaconda
conda activate snakemake
qlogin -l h_vmem=4G
cd ./ne-traj-pipe
snakemake --profile profile > test.log 2>&1 &
```
Otherwise, you can define in your command line. The paramaters is followed by [Configuration](#configuration)
```bash
module load anaconda
conda activate snakemake
qlogin -l h_vmem=4G
cd ./ne-traj-pipe
snakemake --profile profile \
    --config vcffile="test" \
    ref_genome="reference.fa" \
    ancestral_genome="ancestral.fa" \
    workdir="./results" \
    numSampling=5 \
    sampleSize=50 \
    > test.log 2>&1 &
```

## Cluster Usage
set up in `/profile`
### sge (recommend)
```yaml
cluster: "qsub -P roslin_HighlanderLab \
         -o {log_dir}/{rule}.{jobid}.out \
         -e {log_dir}/{rule}.{jobid}.err \
         -l h_vmem={cluster.mem} \
         -l h_rt={cluster.time} \
         -pe sharedmem {cluster.threads} \
         -R y"
```
### drmaa
please check drmaa path in your ```~/.bash_profile``` or equivalent: ```export DRMAA_LIBRARY_PATH=/exports/applications/gridengine/ge-8.6.5/lib/lx-amd64/libdrmaa.so```
Also must login a node, otherwise report error 
```
OpenBLAS blas_thread_init: pthread_create failed for thread 7 of 16: Resource temporarily unavailable
OpenBLAS blas_thread_init: RLIMIT_NPROC 200 current, 200 max
```

check if conda.sh is included in ```~/.bashrc```, if not it will report error
```/bin/bash: conda: command not found```

## Installation
`git clone git@github.com:HighlanderLab/ne-traj-pipe.git`

### Environment Setup
Recommondation of building environments on Eddie
#### PSMC
`source activate ne_pipe_psmc`  
Origianally created as:  
`conda create -n ne_pipe_psmc -c bioconda -c conda-forge psmc`

#### MSMC2 
`source activate ne_pipe_msmc2`  
`msmc2_Linux`
Originally created as:  
`conda create -n ne_pipe_msmc2 -c bioconda -c conda-forge msmc2`

#### SMC++ 
Should be run via singularity
```
module load singularity/3.11.4 # or similar
export SINGULARITY_TMPDIR="/exports/eddie/scratch/username/"
singularity run docker://terhorst/smcpp:latest -h # replace -h by SMC++ parameters
```
test:
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

#### phlash
phlash requires a GPU to run fast. THere are GPUs available on Eddie, but we probably have to play around a little to figure out how to use these. Hannes has set up two versions.

##### No GPU
This is not a conda environment! Python venv in:  
`/exports/cmvm/eddie/eb/groups/HighlanderLab/share/venvs/ne_pipe_phlash_nogpu`
```
##Setup code (so we remember)
#python -m venv --copies ne_pipe_phlash_nogpu
#source ne_pipe_phlash_nogpu/bin/activate
#python -m pip install -U pip setuptools
#python -m pip install 'phlash'
```
##### With GPU support
This is not a conda environment! Python venv in:  
`/exports/cmvm/eddie/eb/groups/HighlanderLab/share/venvs/ne_pipe_phlash_gpu`
```
#Setup code
python -m venv --copies ne_pipe_phlash_gpu
source ne_pipe_phlash_gpu/bin/activate
python -m pip install -U pip setuptools
python -m pip install 'phlash[gpu]'
```

To activate, run:  
```
source /exports/cmvm/eddie/eb/groups/HighlanderLab/share/venvs/ne_pipe_phlash_nogpu/bin/activate
# then, try
python
>>> import phlash
```
#### GONE 
`git clone https://github.com/esrud/GONE.git`

#### Stairway
get easySFS.py 
`git clone https://github.com/isaacovercast/easySFS.git`
Dowload package from https://github.com/xiaoming-liu/stairway-plot-v2
`unzip git_soft/stairway_plot_v2.1.2.zip -d /destination/path/`

#### Relate
Dowload package from https://myersgroup.github.io/relate/
`tar -xvzf ./relate_v1.2.2_x86_64_dynamic.tgz -C /destination/path/`

2. Ensure `conda.sh` is in `~/.bashrc`

#### dadi-cli
`conda create -n ne_pipe_dadiCli -c conda-forge dadi-cli` 

#### fastsimcoal2
`conda create -n ne_pipe_fastsimcoal2 -c bioconda fastsimcoal2`




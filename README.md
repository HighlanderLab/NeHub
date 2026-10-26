# NeHub

A Nextflow pipeline for estimating effective population size (Ne) trajectories from whole-genome sequence data, supporting multiple demographic inference tools.

> **Note:** The NeHub Docker image is currently **private**. It will be made public in a future release. If you wish to use Docker or Singularity, you can either build the image yourself (see [Running with Docker / Singularity](#running-with-docker--singularity)) or contact the author for access to pull the image directly.

> **Chromosome naming:** The pipeline currently only recognises chromosomes in the format `chr1`, `chr2`, ..., `chrX`. Chromosomes named `1`, `2`, etc. (without the `chr` prefix) will not be processed correctly. Please ensure your VCF files and reference genome use the `chr` prefix convention.

> **Docker vs Conda:** If you are not familiar with Docker or Singularity, we recommend using the **Conda-based installation** described in the [Installation](#installation) section. `PSMC` and `MSMC2` are included in `basic_tools.yml`. `SMC++`, `pHlash`, and `tsinfer/tsdate` each require a separate Conda environment as detailed below. `Stairway Plot 2`, `Relate`, and `GONE` are installed as standalone packages, with their installation paths specified in the pipeline configuration.

---

## Table of Contents
- [Available Tools](#available-tools)
- [Configuration](#configuration)
- [Data Processing](#data-processing)
- [Usage](#usage)
- [Cluster Usage](#cluster-usage)
- [Running with Docker / Singularity](#running-with-docker--singularity)
- [Installation](#installation)
- [VCF Header Fix](#vcf-header-fix)
- [Important Notes](#important-notes)

---

## Available Tools

| Tools | Input | Link | Install Approach | Available Now |
|---|---|---|---|---|---|
| GONE | plink file (.ped, .bim) | https://github.com/esrud/GONE | git clone | Yes |
| PSMC | fastq | https://github.com/lh3/psmc | conda install | Yes |
| MSMC2 | vcf | https://github.com/stschiff/msmc2 | conda install | Yes |
| SMC++ | vcf.gz | https://github.com/popgenmethods/smcpp | conda install |Yes |
| pHlash | vcf.gz | https://github.com/jthlab/phlash | pip install | Yes |
| tsinfer and tsdate | ts | https://tskit.dev/ | conda install | Yes |
| Relate | haps and sample | https://myersgroup.github.io/relate/ | package | Yes |
| Stairway Plot v2 | SFS | https://github.com/xiaoming-liu/stairway-plot-v2 | package | Yes |

### Memory Requirements
- SMC++ requires at least 16 GB memory if use `singularity run docker://terhorst/smcpp:latest`
- pHlash benefits from GPU acceleration
- Java memory settings:
```bash
export NXF_OPTS='-Xss256k -Xmx512m'
```

---

## Configuration

The pipeline uses configuration files located in the `conf/` directory, including global parameters and tool-specific configurations.

### Parameter Input Methods

1. **Command Line Parameters**: Direct input through the command line
```bash
nextflow -log ${workdir}/logs/nextflow_${task}_${replicate}.log run ${workdir}/main.nf \
    --chromosomes 1  \
    --population_file popLabel.txt \
    --GONE.installation_path "/path/to/GONE/" \
```

2. **JSON Configuration File**: Partial parameter input through a JSON file

Detail see in `./conf/params.config`. 

```json
{
    "input_dir": "/test/",
    "method": "GONE,PSMC,tskit,SMCpp,Relate,Phlash,stairway,MSMC",
    "chromosomes": [1]
}
```

> **Note:** If use `singularity`, the reference genome, ancestral genome, and container image paths must all be **absolute paths**. Relative paths (e.g. `./tests/test_data/test.fa`) will not work inside Singularity containers as the working directory differs from the host.


## Usage


### Basic Usage 

#### Conda environment
> **Note:** Before starting, please install conda environment, `Stairway Plot 2`, `Relate`, and `GONE`

```bash
module load nextflow
export NXF_OPTS='-Xss256k -Xmx512m'
# or make sure they are iin your environment
nextflow run main.nf -profile local_no_Nehub_docker \
    -params-file params.json \
    --input_dir /path/to/input \
    --output_dir /path/to/output
```

#### Container environment 
> **Note:** Before starting, please pull mirror and build `.sif`. Replace path of `.sif` to container setting of local_singularity in nextflow.config 

```bash
module load singulairty
module load nextflow
export NXF_OPTS='-Xss256k -Xmx512m'
# or make sure they are iin your environment
nextflow run main.nf -profile local_singularity \
    -params-file params.json \
    --input_dir /path/to/input \
    --output_dir /path/to/output
```

### Using a Configuration File

Parameter file refers to`./conf/params.config`, and `./conf/resources.config` has two setting for large resource requirement (`large`) and small resource requirement (`small`). Here, `large` is used for real data, `small` is used for simulation data.

```bash
nextflow run main.nf -profile local_singularity,small \
    -params-file params.json \
    --input_dir /path/to/input \
    --output_dir /path/to/output
```

### Running Multiple Analyses in Parallel
```bash
nextflow run main.nf -profile local_singularity,large \
    -params-file params.json \
    --input_dir /path/to/input \
    --output_dir /path/to/output
```

##### The pipeline provides two resource profiles:

| Profile | Recommended For | Samples | Chromosomes | Memory |
|---|---|---|---|---|
| `small` | Testing, small datasets, simulation | ≤ 100 | ≤ 2 | Reduced (e.g. 4–8 GB per process) |
| `large` | Production runs, large cohorts | > 100 | Many | Full (e.g. 16–160 GB per process) |
| *(none)* | Uses default resources defined inline in each module; resource requirements are similar to `small` | — | — | — |

| Profile | Executor | Software from | Container | Concurrency | Error strategy | Intended use |
|---|---|---|---|---|---|---|---|
| `local_no_Nehub_docker` | local (`queueSize 8`, `maxForks 9`) | Conda envs, activated per label | none | ≤ 8–9 parallel tasks in one job | `ignore` | Eddie or local machine with all conda envs installed |
| `local_singularity` | local (`queueSize 8`, `maxForks 9`) | NeHub Singularity image | same `nehub.sif` (`--nv` GPU flag) | ≤ 8–9 parallel tasks in one job | `ignore` | Eddie, everything runs inside one interactive or array job |
| `sge_singularity` | SGE (`queueSize 8`) | NeHub Singularity image | `/exports/eddie/scratch/$USER/nehub.sif` | 8 queued jobs | `ignore` | Eddie, each process submitted as its own SGE job |
| `sge_no_singularity` | SGE (`queueSize 20`, 10 s submit rate) | Conda envs (Eddie's central anaconda) + `module load` for plink | none (but Singularity is enabled) | 20 queued jobs, `maxForks 9` | `retry` ×3 | Eddie production runs; per-job `h_vmem` = task memory + 4 GB; phlash on the GPU queue |
| `github` | local | Conda envs created in CI | none |  default | default (stop) | GitHub Actions CI; also loads `conf/test.config` |

---

## Cluster Usage

SGE configuration is defined in `nextflow.config`. The pipeline is optimised for SGE cluster environments (e.g. Eddie, University of Edinburgh).

### Submitting on Eddie

```bash
qsub submit_nextflow.sh
```
#### Base on Conda and Installation

```bash
# submit_nextflow.sh
#!/bin/bash
#$ -N nextflow_master
#$ -l h_vmem=8G
#$ -l h_rt=48:00:00
#$ -pe sharedmem 4
#$ -cwd
#$ -j y

module load java
module load nextflow
export NXF_OPTS='-Xms1g -Xmx4g'

nextflow run main.nf \
    -profile local_singularity,large \
    -params-file params.json \
    --input_dir /path/to/input \
    --output_dir /path/to/output \
    --chromosomes 1 \
    --population_file popLabel.txt
```



#### Base on .sif 
```bash
# submit_nextflow.sh
#!/bin/bash
#$ -N nextflow_master
#$ -l h_vmem=8G
#$ -l h_rt=48:00:00
#$ -pe sharedmem 4
#$ -cwd
#$ -j y

module load singularity
module load java
module load nextflow
export NXF_OPTS='-Xms1g -Xmx4g'

export NXF_SINGULARITY_CACHEDIR=/exports/path/to/singularity_cache
mkdir -p $NXF_SINGULARITY_CACHEDIR

nextflow run main.nf \
    -profile local_singularity,large \
    -params-file params.json \
    --input_dir /path/to/input \
    --output_dir /path/to/output \
    --chromosomes 1 \
    --population_file popLabel.txt
```


### Monitor Jobs

```bash
qstat                        # view job status
tail -f .nextflow.log        # view Nextflow log in real time
```

### Resume Interrupted Runs

```bash
nextflow run main.nf -profile eddie_singularity,large -resume
```

---

## Running with Docker / Singularity

> **Note:** If you are not familiar with Docker or Singularity, skip this section and follow the [Conda-based installation](#installation) instead.

Pre-pull the image once to avoid multiple nodes downloading simultaneously.
> **Note:** This step requires at least 16 GB memory. Request an interactive node first:
> ```bash
> qlogin -l h_vmem=16G
> ```

**Option 1: Pull from Docker Hub (requires access credentials)**
```bash
module load singularity
export SINGULARITY_DOCKER_USERNAME=aprilyuzhang
export SINGULARITY_DOCKER_PASSWORD=${PASSWORD} # contact author if need
export SINGULARITY_TMPDIR=/exports/eddie/scratch/$USER/tmp
mkdir -p /exports/eddie/scratch/$USER/tmp

singularity pull \
    --tmpdir /exports/eddie/scratch/$USER/tmp \
    ./nehub.sif \
    docker://aprilyuzhang/nehub:latest
```
**Option 2: Build from Dockerfile (no credentials required)**
```bash
# 1. Build Docker image locally on your machine
docker build --build-arg BUILDPLATFORM=linux/amd64 \
    -t nehub:latest \
    -f docker/Dockerfile .

# 2. Save as tar
docker save nehub:latest -o nehub.tar
gzip nehub.tar

# 3. Transfer to Eddie
scp nehub.tar.gz /path/to/image

# 4. On Eddie, convert to SIF
module load singularity
# set SINGULARITY_TMPDIR to avoid FATAL ERROR: Out of memory (cache_alloc)
export SINGULARITY_TMPDIR=/exports/eddie/scratch/$USER/tmp
mkdir -p $SINGULARITY_TMPDIR
gunzip /path/to/image/nehub.tar.gz
singularity build \
    ./nehub.sif \
    docker-archive://path/to/image/nehub.tar
```

**Option 3: Build from Dockerfile and push to private repository (no credentials required)**
```bash
# 1. Build Docker image locally on your machine
docker build --build-arg BUILDPLATFORM=linux/amd64 \
    -t nehub:latest \
    -f docker/Dockerfile .

docker tag nehub:latest aprilyuzhang/nehub:latest
docker push aprilyuzhang/nehub:latest
module load singularity
export SINGULARITY_DOCKER_USERNAME=aprilyuzhang
export SINGULARITY_DOCKER_PASSWORD=${PASSWORD}
export SINGULARITY_TMPDIR=/exports/eddie/scratch/$USER/tmp
mkdir -p /exports/eddie/scratch/$USER/tmp

singularity pull \
    --tmpdir /exports/eddie/scratch/$USER/tmp \
    $USER/nehub.sif \
    docker://aprilyuzhang/nehub:latest
```

---

## Installation

> **Recommendation:** If you are not familiar with Docker, use the Conda-based installation below. Each tool is installed into a separate Conda environment to avoid dependency conflicts.

### Requirements

- Nextflow >= 24.10.2
- Conda (recommended: [Miniconda](https://docs.conda.io/en/latest/miniconda.html))
- Java >= 17

### Clone the Repository

```bash
git clone https://github.com/AprilYUZhang/NeTrajHub.git
cd NeTrajHub
```

### Conda Environments
| Conda env | YAML file | Main contents | Used by process label | Notes |
|---|---|---|---|---|
| `basic_Ne` | `conf/envs/basic_tools.yml` | PSMC, MSMC2, plink, bcftools, other basic tools | `basic_Ne`, `plink` | Core environment for most steps |
| `ts_tools` | `conf/envs/ts_tools.yml` | tskit, tsinfer, tsdate, tszip, msprime | `ts_tools` | Then run `conda install -n ts_tools -y -c conda-forge tszip msprime` |
| `bcftools` | `conf/envs/bcftools.yml` | bcftools | `bcftools` | VCF handling and consensus FASTA |
| `smcpp` | `conf/envs/smcpp.yml` | SMC++ | `smcpp` | Requires Python 3.10 (`longintrepr.h` was removed in Python 3.11+) |
| `phlash` | `conf/envs/phlash.yml` | phlash | `ne_pipe_phlash` | Works with both GPU and CPU builds |

> **Note:** if smcpp Installation has issue, alternatively, use the Singularity image (no installation required and more stable)

> **Note:** Ensure your VCF file header includes chromosome lengths before running SMC++. See [VCF Header Fix](#vcf-header-fix) below.

> **Note:** If no GPU is detected, phlash automatically falls back to CPU. Installing the standard version (`pip install phlash`) will use GPU if available, or CPU otherwise — no need to maintain two separate environments.

##### Known potential pHlash Installation Issues

**1. `dinopy` build failure:**

If you encounter a Cython compilation error related to `isinstance(seq, (int, long))`:
```bash
# download dinopy source
pip download dinopy --no-deps -d ./dinopy_src
# unpack and edit processors.pyx:
# change: isinstance(seq, (int, long))
# to:     isinstance(seq, int)
pip install ./dinopy-3.0.0/
```

**2. JAX >= 0.6.0 incompatibility:**

If JAX version is >= 0.6.0, edit `gpu.py` in the phlash installation:
```python
# change from:
return jax.pure_callback(
    kern, result_shape_dtype, pp=params, index=index, grad=grad,
    vmap_method='legacy_vectorized'
)

# to:
return jax.pure_callback(
    kern, result_shape_dtype, pp=params, index=index, grad=grad,
)
```

#### GONE

```bash
git clone https://github.com/esrud/GONE.git /path/to/GONE
chmod -R +x /path/to/GONE/Linux/PROGRAMMES/
```

#### Stairway Plot v2

```bash
git clone https://github.com/isaacovercast/easySFS.git /path/to/easySFS
chmod +x /path/to/easySFS/easySFS.py

# Download stairway_plot_v2.1.3.zip from:
# https://github.com/xiaoming-liu/stairway-plot-v2
unzip stairway_plot_v2.1.3.zip
```

#### Relate

Download the pre-built binary from https://myersgroup.github.io/relate/:
```bash
tar -xvzf relate_v1.2.2_x86_64_dynamic.tgz -C /path/to/relate/
```

Or build from source:
```bash
git clone --depth 1 https://github.com/MyersGroup/relate.git /path/to/relate_src
cd /path/to/relate_src
rm -rf build && mkdir build && cd build
cmake .. && make -j$(nproc)
cp /path/to/relate_src/bin/* /path/to/relate/
```

---

## VCF Header Fix

Some tools (e.g. SMC++) require the VCF header to include chromosome lengths. Use the following script to add them:

```bash
i=chr_number
bcftools view -h chr${i}.vcf.gz > old_header.txt

awk 'BEGIN{while(getline < "ref.fa.fai"){len["chr"$1]=$2}}
/^##contig=<ID=/ {match($0,/ID=([^,>]+)/,a); chr=a[1];
    if(chr in len){print "##contig=<ID="chr",length="len[chr]">"}
    else {print $0}}
!/##contig=<ID=/ {print $0}' old_header.txt > new_header.txt

bcftools reheader -h new_header.txt -o chr${i}_new.vcf.gz chr${i}.vcf.gz

mv chr${i}.vcf.gz chr${i}_old.vcf.gz
mv chr${i}_new.vcf.gz chr${i}.vcf.gz
```

---

## Important Notes

- **Chromosome naming:** The pipeline currently only recognises chromosomes in the format `chr1`, `chr2`, ..., `chrX`. Please ensure all input VCF files, reference genomes, and ancestral genomes use the `chr` prefix. Chromosomes named `1`, `2`, etc. will cause the pipeline to fail silently.
- **Docker image:** The Docker image `aprilyuzhang/nehub` is currently private. Conda-based installation is recommended for users without Docker Hub access.
- **GPU support:** pHlash benefits significantly from GPU acceleration. On Eddie, request a GPU node with `-q gpu -l gpu=1` for pHlash jobs.
- **Refernece genome:** Make sure `Refernece genome` is consistent with vcf REF
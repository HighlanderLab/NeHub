#!/usr/bin/env python
"""
Infer a tree sequence from a phased VCF file using tsinfer.

Steps:
  1. Read the VCF and ancestral FASTA.
  2. Build a tsinfer SampleData object (one population per population_group).
  3. Run tsinfer to produce a tree sequence.
  4. Simplify and compress the result to <output>.trees.tsz.

Usage:
    python ts_prepare.py --vcf <vcf.gz> --ancestral <ref.fa> \
        --popLabel <poplabel.txt> --output <prefix> \
        [--max_file_size N] [--chunk_size N] [--thread N]

Note: the VCF must be phased; each sample must appear in the pop-label file.
"""

import cyvcf2
from tqdm import tqdm
import tskit
import tsinfer
import os
from Bio import SeqIO
import argparse
import gzip
import tszip

#https://tskit.dev/tsinfer/docs/stable/tutorial.html#reading-a-vcf
def add_diploid_sites(vcf, samples,ancestral_fa,sample_name):
    """
    Read the sites in the vcf and add them to the samples object.
    """
    # You may want to change the following line, e.g. here we allow
    # "*" (a spanning deletion) to be a valid allele state
    allele_chars = set("ATGCatgc*")
    pos = 0
    progressbar = tqdm(total=samples.sequence_length, desc="Read VCF", unit='bp')
    non_match_AA = 0
    # Get all sample names
    all_samples = vcf.samples

    # Get indices for samples of interest
    sample_indices = []
    for sample in sample_name:
        if sample in all_samples:
            sample_indices.append(all_samples.index(sample))
        else:
            print(f"Error: Sample '{sample}' not found in VCF")
            exit()

    for variant in vcf:  # Loop over variants, each assumed at a unique site
        anc_record=ancestral_fa[str(variant.CHROM)]
        ancestral = anc_record.seq[variant.POS-1]
        progressbar.update(variant.POS - pos)
        if pos == variant.POS:
            print(f"Duplicate entries at position {pos}, ignoring all but the first")
            continue
        else:
            pos = variant.POS
        if any([not phased for _, _, phased in variant.genotypes]):
            raise ValueError("Unphased genotypes for variant at position", pos)
        alleles = [variant.REF.upper()] + [v.upper() for v in variant.ALT]
        # ancestral = variant.INFO.get("AA", ".")  # "." means unknown
        # some VCFs (e.g. from 1000G) have many values in the AA field: take the 1st
        # ancestral = ancestral.split("|")[0].upper()

        if ancestral == "." or ancestral == "":
            ancestral_allele = tskit.MISSING_DATA

        elif ancestral not in alleles:
            alleles.append(ancestral)
            ancestral_allele = alleles.index(ancestral)
            non_match_AA+=1
        # alternatively, you could specify `ancestral = variant.REF.upper()`
        else:
            ancestral_allele = alleles.index(ancestral)
        # Check we have ATCG alleles
        for a in alleles:
            if len(set(a) - allele_chars) > 0:
                print(f"Ignoring site at pos {pos}: allele {a} not in {allele_chars}")
                continue
        # Map original allele indexes to their indexes in the new alleles list.
        genotypes = []
        for sample_index in sample_indices:
            g = variant.genotypes[sample_index]
            genotypes.extend([g[0],g[1]])
        samples.add_site(pos, genotypes, alleles, ancestral_allele=ancestral_allele)

    print("The allele of VCF records excluded for 'AA':", non_match_AA)

def add_pop_indi(popLabel,samples):
    count=0
    for pop in popLabel:
        samples.add_population(metadata={"name": pop})
    sample_name=[]
    for pop in popLabel:
        for i in popLabel[pop]:
            samples.add_individual(ploidy=2, population=count, metadata={"name": i})
            sample_name.append(i)
        count+=1
    return sample_name

def chromosome_ID(vcf):
    assert len(vcf.seqnames) == 1
    return vcf.seqnames[0]

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description="Process some integers.")

# Add arguments
parser.add_argument('--vcf', type=str, required=True, help='Path to the VCF file')
parser.add_argument('--ancestral', type=str, required=True, help='Path to the ancestral genome fa file')
parser.add_argument('--max_file_size', type=int, default=4*500*1024*1024 ,help='max_file_size of Zarr')
parser.add_argument('--chunk_size', type=int, default = 2*1024, help='chunk_size of Zarr')
parser.add_argument('--output', type=str, required=True, help='Path to the output file')
parser.add_argument('--popLabel', type=str, required=True, help='Path to the popLabel file')
parser.add_argument('--thread', type=int, default= 0 , help='Set the threads for tsinfer')
# Parse arguments
args = parser.parse_args()

input = args.vcf

# build new folder for samples
directory_path=os.path.dirname(input)
basename=os.path.basename(input)
basename = basename.removesuffix(".vcf.gz")

vcf_location = input # "ccon0.vcf.gz"
    # NB: could also read from an online version by setting vcf_location to
    # "https://github.com/tskit-dev/tsinfer/raw/main/docs/_static/P_dom_chr24_phased.vcf.gz"

vcf = cyvcf2.VCF(vcf_location)

# ancestral genome
fasta_file = args.ancestral # "reference_genome.fa"
# Read the FASTA file
ancestral_fa={}
if fasta_file.endswith(".gz"):
    with gzip.open(fasta_file, "rt") as handle: # <--- This is crucial
        for record in SeqIO.parse(handle, "fasta"):
            ancestral_fa[record.id]=record
else:
    for record in SeqIO.parse(fasta_file, "fasta"):
        ancestral_fa[record.id]=record

# population label
popLabel= args.popLabel # "ccon0.poplabels"
poplabel_dic={}
with open(popLabel, "r") as f:
    f.readline()
    line=f.readline()
    while line !="":
        li=line.strip("\n").split("\t")
        indi=li[0]
        pop=li[1]
        group=li[2]
        pop_group=pop+"_"+group
        if pop_group not in poplabel_dic:
            poplabel_dic[pop_group]=[indi]
        else:
            poplabel_dic[pop_group].append(indi)
        line=f.readline()

    # "ccon0.samples"

with tsinfer.SampleData(
        path= f"{args.output}.samples", sequence_length=len(ancestral_fa[chromosome_ID(vcf)].seq),max_file_size=args.max_file_size,chunk_size=args.chunk_size,
        ) as samples:
        sample_name = add_pop_indi(poplabel_dic,samples)
        add_diploid_sites(vcf, samples, ancestral_fa,sample_name)
    
print(
        "Sample file created for {} samples ".format(samples.num_samples)
        + "({} individuals) ".format(samples.num_individuals)
        + "with {} variable sites.".format(samples.num_sites),
        flush=True,
    )

ts = tsinfer.infer(samples,num_threads=args.thread)

print("Inferred tree sequence: {} trees over {} Mb ({} edges)".format(
        ts.num_trees, ts.sequence_length / 1e6, ts.num_edges
                ))
ts=ts.simplify()
print("simplified tree sequence: {} trees over {} Mb ({} edges)".format(
        ts.num_trees, ts.sequence_length / 1e6, ts.num_edges
                ))
tszip.compress(ts,f"{args.output}.trees.tsz")
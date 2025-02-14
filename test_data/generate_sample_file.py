import cyvcf2
from tqdm import tqdm
import tskit
from Bio import SeqIO

#https://tskit.dev/tsinfer/docs/stable/tutorial.html#reading-a-vcf
def add_diploid_sites(vcf, samples,ancestral_fa):
    """
    Read the sites in the vcf and add them to the samples object.
    """
    # You may want to change the following line, e.g. here we allow
    # "*" (a spanning deletion) to be a valid allele state
    allele_chars = set("ATGCatgc*")
    pos = 0
    progressbar = tqdm(total=samples.sequence_length, desc="Read VCF", unit='bp')
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
            # alternatively, you could specify `ancestral = variant.REF.upper()`
        else:
            ancestral_allele = alleles.index(ancestral)
        # Check we have ATCG alleles
        for a in alleles:
            if len(set(a) - allele_chars) > 0:
                print(f"Ignoring site at pos {pos}: allele {a} not in {allele_chars}")
                continue
        # Map original allele indexes to their indexes in the new alleles list.
        genotypes = [g for row in variant.genotypes for g in row[0:2]]
        samples.add_site(pos, genotypes, alleles, ancestral_allele=ancestral_allele)

def add_pop_indi(popLabel,samples):
    count=0
    for pop in popLabel:
        samples.add_population(metadata={"name": pop})
        for i in popLabel[pop]:
            samples.add_individual(ploidy=2, population=count, metadata={"name": i})
        count+=1

def chromosome_length(vcf):
    assert len(vcf.seqlens) == 1
    return vcf.seqlens[0]

import argparse

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description="Process some integers.")

# Add arguments
parser.add_argument('--vcf', type=str, required=True, help='Path to the VCF file')
parser.add_argument('--ancestal', type=str, required=True, help='Path to the ancestral genome fa file')
parser.add_argument('--popLabel', type=str, required=True, help='Path to the popLabel file')
parser.add_argument('--output', type=str, required=True, help='Path to the output file')

# Parse arguments
args = parser.parse_args()


vcf_location = args.vcf # "ccon0.vcf.gz"
# NB: could also read from an online version by setting vcf_location to
# "https://github.com/tskit-dev/tsinfer/raw/main/docs/_static/P_dom_chr24_phased.vcf.gz"

vcf = cyvcf2.VCF(vcf_location)

# ancestral genome
fasta_file = args.ancestal # "reference_genome.fa"
SeqIO.parse(fasta_file, "fasta")
# Read the FASTA file
ancestral_fa={}
for record in SeqIO.parse(fasta_file, "fasta"):
    ancestral_fa[record.id]=record

# population label
popLabel= args.popLabel # "ccon0.popLabel"
poplabel_dic={}
with open(popLabel, "r") as f:
    f.readline()
    line=f.readline()
    while line !="":
        li=line.strip("\n").split("\t")
        indi=li[0]
        pop=li[1]
        if pop not in poplabel_dic:
            poplabel_dic[pop]=[indi]
        else:
            poplabel_dic[pop].append(indi)
        line=f.readline()

# "ccon0.samples"

with tsinfer.SampleData(
    path= args.output, sequence_length=chromosome_length(vcf),max_file_size=4*500*1024*1024,chunk_size=2*1024,
) as samples:
    add_pop_indi(poplabel_dic,samples)
    add_diploid_sites(vcf, samples, ancestral_fa)
    
print(
    "Sample file created for {} samples ".format(samples.num_samples)
    + "({} individuals) ".format(samples.num_individuals)
    + "with {} variable sites.".format(samples.num_sites),
    flush=True,
)

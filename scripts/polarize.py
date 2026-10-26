#!/usr/bin/env python
"""
Polarize a VCF file so that REF alleles match the ancestral genome.

For each variant:
  - If REF == ancestral allele: leave unchanged, set INFO/AA=REF.
  - If ancestral allele is in ALT: swap REF and the matching ALT, recode
    genotypes accordingly, set INFO/AA=ancestral.
  - If ancestral allele is absent from REF/ALT: mark FILTER=MISSING,
    set AA='.'.

The output VCF (gzip-compressed) can be filtered with bcftools
  bcftools view -e 'FILTER="MISSING"'
to remove sites with unknown ancestry before SFS computation.

Usage:
    python polarize.py -v <input.vcf[.gz]> -a <ancestral.fa[.gz]> -o <output.vcf.gz>
"""

from Bio import SeqIO
import gzip
import argparse

args = argparse.ArgumentParser(description="Polarize a VCF file based on an ancestral genome.")
args.add_argument("-v", "--vcf", required=True, help="Input VCF file")
args.add_argument("-a", "--anc", required=True, help="Ancestral genome in FASTA format")
args.add_argument("-o", "--out", required=True, help="Output VCF file")
args = args.parse_args()

fasta_file = args.anc
ancestral_fa={}
if fasta_file.endswith(".gz"):
    with gzip.open(fasta_file, "rt") as handle: # <--- This is crucial
        for record in SeqIO.parse(handle, "fasta"):
            ancestral_fa[record.id]=record
else:
    for record in SeqIO.parse(fasta_file, "fasta"):
        ancestral_fa[record.id]=record
vcffile = args.vcf
outfile = args.out

## use gzip? 
read_mode = 'r'
if vcffile.endswith(".gz"):
    read_mode = 'rt'
    ofunc = gzip.open
else:  
    ofunc = open

infile = ofunc(vcffile, read_mode)
with ofunc(vcffile, read_mode) as fin, gzip.open(outfile, 'wt') as fout:
    has_AA_info = False
    last_line = ''
    for line in fin:
        if line.startswith('##INFO=<ID=AA'):
            has_AA_info = True

        if last_line.startswith('##INFO') and not line.startswith('##INFO'):
            # Add the AA info line if it doesn't exist
            if not has_AA_info:
                fout.write('##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele">' + '\n')

        if last_line.startswith('##') and not line.startswith('##'):
            fout.write('##FILTER=<ID=MISSING,Description="Contains missing ancestral allele">'+ '\n')

        if line.startswith('#'):
            fout.write(line)
            last_line = line
            continue
        
        cols = line.rstrip('\n').split('\t')
        info_field = cols[7]
        ref = cols[3]
        alt = cols[4]
        chr = cols[0]
        pos = cols[1]
        anc = ancestral_fa[chr].seq[int(pos) - 1]

        if ref != anc:
            # Polarize the REF and ALT alleles
            alt_li = alt.split(",")
            if anc in alt_li:
                # Find which ALT allele matches ancestral
                anc_alt_index = alt_li.index(anc) + 1  # +1 because ALT indices start at 1

                # Swap REF and ALT
                cols[3] = anc  # Set ancestral as new REF
                alt_li[anc_alt_index - 1] = ref  # Put old REF in ALT position
                cols[4] = ",".join(alt_li)

                for i in range(9, len(cols)):
                    genotype = cols[i]
                    if genotype == './.' or genotype == '.|.':
                        continue
                    separator = '|' if '|' in genotype else '/'
                    alleles = genotype.split(separator)

                    new_alleles = []
                    for allele in alleles:
                        if allele == '.':
                            new_alleles.append('.')
                        elif allele == '0':  # REF -> ALT position
                            new_alleles.append(str(anc_alt_index))
                        elif allele == str(anc_alt_index):  # ALT -> REF (0)
                            new_alleles.append('0')
                        else:
                            new_alleles.append(allele)  # Other ALTs unchanged
                    cols[i] = separator.join(new_alleles)

            else:
                cols[6] = "MISSING"
                anc = '.'


        # Find AA=... in INFO and replace it with AA=<REF>
        if 'AA=' in info_field:
            # Extract current AA value
            info_items = info_field.split(";")
            new_items = []
            for item in info_items:
                if item.startswith('AA='):
                    new_items.append(f"AA={anc}")
                else:
                    new_items.append(item)
            info_field = ";".join(new_items)
        else:
            info_field = info_field + ";AA=" + anc
        cols[7] = info_field

        fout.write('\t'.join(cols) + '\n')

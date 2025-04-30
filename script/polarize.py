from Bio import SeqIO
import gzip
import argparse
args= argparse.ArgumentParser(description="Polarize a VCF file based on a reference genome.")
args.add_argument("-v", "--vcf", required=True, help="Input VCF file")
args.add_argument("-a", "--anc", required=True, help="Ancestral genome in FASTA format")
args.add_argument("-o", "--out", required=True, help="Output VCF file")
args = args.parse_args()
fa=SeqIO.parse(args.anc, "fasta")
ancestral_fa={}
for record in fa:
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
                cols[3] = anc
                alt_index = alt_li.index(anc)
                alt_li[alt_index] = ref
                cols[4] = ",".join(alt_li)
                genotypes = cols[9:]
                for genotype in genotypes:
                    if genotype != '.':
                        gt = genotype.split(":")
                        gt[0] = gt[0].replace(alt_index+1, "M")
                        gt[0] = gt[0].replace("0", alt_index+1)
                        gt[0] = gt[0].replace("M", "0")
                        genotype = ":".join(gt)

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

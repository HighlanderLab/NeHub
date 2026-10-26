#!/usr/bin/env python
"""
Generate minimal synthetic test data for all 8 Ne-estimation methods.

Outputs (written to --outdir):
  test_chr1.vcf.gz          phased VCF, 8 diploid individuals (2 pop × 4),
                             bgzip-compressed and tabix-indexed
  test_chr1.vcf.gz.tbi      tabix index
  test.fa                   ancestral/reference FASTA (chr id = "chr1", 100 kb)
  popLabel.txt              population label file (TSV: sample_id/population/group/sex)
  test.chr1.genetic_map.txt flat-rate genetic map for Relate

Method coverage:
  GONE     → VCF → plink (pipeline handles conversion)
  SMCpp    → VCF → smc format (pipeline handles conversion)
  stairway → VCF + ancestral genome (polarization) + popLabel
  tskit    → VCF + ancestral genome + popLabel
  PSMC     → VCF + reference genome (consensus calling)
  Relate   → VCF + ancestral genome + genetic map + popLabel
  Phlash   → VCF + popLabel
  MSMC     → VCF + popLabel (pipeline handles tab conversion)

Dependencies (pip install):
  msprime numpy pandas
System tools:
  bgzip, tabix  (from htslib; install with: sudo apt-get install tabix)
# msprime-1.2.0_tskit-0.5.2_pyslim-1.0
"""

import argparse
import subprocess
import io
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import msprime
import random


# ─────────────────────────────────────────────────────────────────────────────

def simulate(n_per_pop: int, seq_length: int, seed: int):
    """Simulate a two-population tree sequence using msprime."""
    demography = msprime.Demography()
    demography.add_population(name="pop0", initial_size=10_000)
    demography.add_population(name="pop1", initial_size=10_000)
    demography.add_population(name="ancestral", initial_size=20_000)
    demography.add_population_split(
        time=500, derived=["pop0", "pop1"], ancestral="ancestral"
    )

    ts = msprime.sim_ancestry(
        samples={"pop0": n_per_pop, "pop1": n_per_pop},
        demography=demography,
        recombination_rate=1e-8,
        sequence_length=seq_length,
        random_seed=seed,
    )
    return msprime.sim_mutations(ts, rate=1e-8, random_seed=seed)


def write_vcf(mts, outdir: Path):
    """Write bgzipped, tabix-indexed VCF. File: test_chr1.vcf.gz"""
    raw_vcf = outdir / "test_chr1.vcf"
    buf = io.StringIO()
    mts.write_vcf(buf, contig_id="chr1",position_transform=lambda x: np.array(x) + 1)
    
    content = buf.getvalue()
    
    # check content is not empty
    if not content:
        raise ValueError("VCF content is empty — check mts has mutations")
    print(f"  VCF content: {len(content)} bytes, {mts.num_mutations} mutations")
    
    raw_vcf.write_text(content)

    # verify file was written
    if raw_vcf.stat().st_size == 0:
        raise ValueError(f"VCF file is empty: {raw_vcf}")

    gz_path = str(raw_vcf) + ".gz"
    # module load igmm/apps/tabix/0.2.6 
    subprocess.run(["bgzip", "-f", str(raw_vcf)], check=True)
    subprocess.run(["tabix", "-p", "vcf", gz_path], check=True)
    print(f"  VCF : {gz_path} ({mts.num_mutations} variants)")

def write_fasta(ts,seed: int, outdir: Path):
    """Write a random ancestral/reference FASTA.  Chromosome id = 'chr1'.

    Both the tskit, Relate, PSMC, and stairway modules need a FASTA whose
    sequence id matches the CHROM column in the VCF ('chr1').
    """
    print(f"  DEBUG outdir type: {type(outdir)}, value: {outdir}")
    print(f"  DEBUG ts type: {type(ts)}")
    print(f"  DEBUG sequence_length: {ts.sequence_length}")
    print(f"  DEBUG num_sites: {ts.num_sites}")
    print(f"  DEBUG num_mutations: {ts.num_mutations}")
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    try:
        s = [random.choices(
                ["A", "T", "C", "G"],
                [0.225, 0.225, 0.275, 0.275],
                k=1)[0]
             for i in range(int(ts.sequence_length))]
        print(f"  DEBUG sequence generated: {len(s)} bases")
    except Exception as e:
        print(f"  ERROR generating sequence: {e}")
        raise
    for i in ts.sites():
        s[int(i.position)]=i.ancestral_state
    # write reference genome, non-variant sites will be written randomly (ATCG)
    out = outdir / "test.fa"
    f=open(out,"w")
    f.write(f">chr1\n")
    for i in range(len(s)):
        if (i%60==0) and (i!=0):
            f.write("\n")
        f.write(s[i])
    f.write("\n")
    f.close()
    print(f"  FASTA: {outdir}/test.fa")

def write_poplabel(n_per_pop: int, outdir: Path):
    """Write TSV population label.  msprime VCF sample names = tsk_0 … tsk_N."""
    rows = [
        (f"tsk_{i}", "pop0", "pop0", "NA") for i in range(n_per_pop)
    ] + [
        (f"tsk_{i}", "pop1", "pop1", "NA")
        for i in range(n_per_pop, 2 * n_per_pop)
    ]
    df = pd.DataFrame(rows, columns=["sample_id", "population", "group", "sex"])
    out = outdir / "popLabel.txt"
    df.to_csv(out, sep="\t", index=False)
    print(f"  popLabel: {out} ({len(df)} samples)")


def write_genetic_map(seq_length: int, outdir: Path):
    """Write a flat-rate genetic map for Relate.

    Format expected by Relate --mode MakeChunks:
        pos  COMBINED_rate  Genetic_Map(cM)

    File name must match the pattern used in Relate1.nf:
        <prefix>.chr<N>.genetic_map.txt
    """
    positions = list(range(0, seq_length + 1, 1000))
    cM_per_bp = 1e-6        # 1 cM / Mb
    out = outdir / "test.chr1.genetic_map.txt"
    cumulative = 0.0
    prev = 0
    lines = ["pos COMBINED_rate Genetic_Map(cM)\n"]
    for pos in positions:
        if pos > 0:
            cumulative += (pos - prev) * cM_per_bp
        lines.append(f"{pos} 1.0 {cumulative:.8f}\n")
        prev = pos
    out.write_text("".join(lines))
    print(f"  Genetic map: {out}")


# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outdir",    default="tests/test_data",
                        help="Output directory (created if absent)")
    parser.add_argument("--n_per_pop", type=int, default=4,
                        help="Diploid individuals per population (default: 4)")
    parser.add_argument("--seq_length", type=int, default=1_000_000,
                        help="Simulated sequence length in bp (default: 100 000)")
    parser.add_argument("--seed",      type=int, default=42,
                        help="Random seed (default: 42)")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Simulating {2 * args.n_per_pop} individuals "
          f"({args.n_per_pop} per population), "
          f"{args.seq_length:,} bp, seed={args.seed}")

    mts = simulate(args.n_per_pop, args.seq_length, args.seed)

    write_vcf(mts, outdir)
    write_fasta(mts, args.seed, outdir)
    write_poplabel(args.n_per_pop, outdir)
    write_genetic_map(args.seq_length, outdir)

    print(f"\nTest data ready in: {outdir.resolve()}")


if __name__ == "__main__":
    main()

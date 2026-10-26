#!/usr/bin/env python
"""
Date a tree sequence produced by ts_prepare.py using tsdate.

Reads a compressed tree sequence (<input>.trees.tsz), runs tsdate with
the given mutation rate, and writes the dated tree sequence to
<input>.trees.dated.tsz.

Usage:
    python ts_date.py --input <prefix>.trees.tsz --mutation_rate <rate>
"""

import tsdate
import argparse
import tszip

parser = argparse.ArgumentParser(description="Date a tree sequence with tsdate.")

# Add arguments
parser.add_argument('--input', type=str, required=True, help='Path to the sampledata file')
parser.add_argument('--mutation_rate', type=float, required=True, help='mutation rate')
# parser.add_argument('--output', type=str, required=True, help='Path to the output trees file')
# Parse arguments
args = parser.parse_args()

input_tree=args.input
mutation = args.mutation_rate
output = input_tree.replace(".trees.tsz",".trees.dated.tsz")
ts=tszip.load(input_tree)
dated_ts = tsdate.date(ts, mutation_rate=mutation)
tszip.compress(dated_ts,output)

import re
import csv

import argparse

from pathlib import Path
import pandas as pd

from Bio import SeqIO
from classes import Sequence, FastxFile

parser = argparse.ArgumentParser(description='''Script to extract read belonging to a construct in constrcuts.csv to a fastq file for single plasmid assembly.''', 
	prog='python3 parts2fasta.py')
parser.add_argument('reads', metavar='IN.fastq.gz', type=Path,
	help='Raw reads.')
parser.add_argument('constructs', metavar='constructs.csv', type=Path,
	help='Constructs file.')
parser.add_argument('output', metavar='OUT', type=Path,
	help='Output folder.')

parser.add_argument('-p', '--percent', default=0.01, type=float, help='bar help')
parser.add_argument('-c', '--count', default=0, type=int, help='bar help')

args = parser.parse_args()

reads_group = {}
constructs = pd.read_csv(args.constructs)
for read, group in constructs.groupby('read'):
	name = [f'{g.TU}_{g.name}' for g in group.itertuples()]

	reads_group.setdefault('_'.join(name), []).append(read)

cutoffN = max(int(constructs.shape[0] * args.percent), args.count)
print(f'Cutoff at: {cutoffN}')

groups = {}
files = {}

for group, reads in reads_group.items():
	if len(reads) > cutoffN:
		print(f'BC{len(files) + 1}: {group} - {len(reads)} reads')
		barcode = f"barcode{len(files):02.0f}"

		(args.output / barcode).mkdir(exist_ok=True)
		files[group] = (barcode, len(reads), [], open(args.output / barcode / 'reads0.fastq', 'w'))

		for read in reads:
			read = read.split('_')[0]
			groups[read] = group

with FastxFile(args.reads) as f:
	reads = SeqIO.parse(f, 'fastq')

	for read in reads:
		if read.id in groups:
			meta = files[groups[read.id]]
			meta[2].append(len(read.seq))
			SeqIO.write(read, meta[3], "fastq")

with open(args.output / "sample_sheet.csv", 'w') as f:
	f.write("alias,barcode,type,approx_size\n")

	for group, meta in files.items():
		meta[3].close()
		sample = group.replace('|', '_')
		f.write(f"{sample},{meta[0]},test_sample,{int(sum(meta[2])/meta[1])}\n")
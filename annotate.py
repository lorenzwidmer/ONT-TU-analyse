import subprocess
import argparse
import tempfile
import shutil

from pathlib import Path
from collections import namedtuple
from zipfile import ZipFile
from io import TextIOWrapper

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

import pandas as pd
import numpy as np

Sequence = namedtuple('Sequence', ['name', 'seq', 'alignments'])
Alignment = namedtuple('Alignment', ['name', 'part', 'interval', 'coverage', 'score'])

parser = argparse.ArgumentParser(description='''Script to find muCAR receptors in NanoPore Sequencing files.''', 
	prog='python3 -m sbo-builder')
parser.add_argument('fastqs', nargs='+', metavar='seq.fastq', type=Path,
	help='File containing the raw sequencing data.')
parser.add_argument('-o', '--output', type=Path, default=Path('./'),
	help='Name of the output file.')
parser.add_argument('-a', '--annotate', dest='annotate', action='store_true',
	help='Create zip file with all annotated sequences.')
parser.add_argument('-p', '--parts', dest='parts', type=str, default='start,promo,rbs,signal,ntag(o),cds,ctag,stop,term,end(o)',
	help='Parts required for a valid receptor.')
parser.add_argument('-p', '--primer', dest='parts', type=str, default='gaattcgcggccgcttctagag,AAAAA',
	help='Parts required for a valid receptor.')

args = parser.parse_args()

def overlap(a, b):
	return min(a.right, b.right) - max(a.left, b.left) + 1

def makeFeatureLocation(alignment):
	if alignment.interval.left < 0:
		return FeatureLocation(-alignment.interval.right, -alignment.interval.left, strand=-1)
	return FeatureLocation(alignment.interval.left, alignment.interval.right, strand=1)

def reorderAlignments(alignments, order):
	if len(alignments) < 1:
		return []
	if len(alignments) > 1 and alignments[0].name == alignments[-1].name and alignments[0].part == alignments[-1].part:
		alignments.pop()

	parts = [a.part.split('-')[0] for a in alignments]
	start = np.argsort(np.array([order.get(p, 100) for p in parts]))[0]
	
	return alignments[start:] + alignments[:start]

root = Path(__file__).parent
folder = root / 'blast'

headers = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qcovs']
primerF, primerR = args.primer.split(',')

sequences = {}
for fastq in args.fastqs:
	with fastq.open('r') as f:
		fasta_sequences = SeqIO.parse(f, 'fastq')

		for fasta in fasta_sequences:
			start = max(0, fasta.find('primerF'))
			end = max(0, fasta.find('primerR'))

			if start == 0 and end == 0:
				print(f'Skipped {fasta.id}')
				continue
			elif end < start:
				newy = fasta[start:] + fasta[0:end]
			else:
				newy = fasta[start:end]

			fasta.id = f'{fasta.id}_{len(sequences)}'
			fasta.annotations["molecule_type"] = "DNA"
			sequences[fasta.id] = Sequence(fasta.id, fasta, [])

with tempfile.TemporaryDirectory() as tempdir:
	tempdir = Path(tempdir)
	with (tempdir / 'plasmids.fasta').open("w") as f:
		SeqIO.write([sequence.seq for sequence in sequences.values()], f, "fasta")

	if not (tempdir / 'plasmids.ndb').is_file():
		subprocess.run([
			'makeblastdb', '-in', 'plasmids.fasta', '-dbtype', 'nucl', '-parse_seqids', '-out', 'plasmids', '-title', 'no-idea'
		], cwd=tempdir)

		subprocess.run([
			'blastn', '-out', tempdir / 'alignments.csv', '-outfmt', ' '.join(['6'] + headers),
			'-max_target_seqs', str(len(sequences)),
			'-query', root / 'parts.fasta', '-db', tempdir / 'plasmids'
		], cwd=tempdir)

	results = pd.read_csv(tempdir / 'alignments.csv', sep="\t", names=headers)
	shutil.copy(tempdir / 'alignments.csv', args.output / 'alignments.csv')

for i, row in results.iterrows():
	name, part = row['qseqid'].split('|')
	modi = -2 * (row['sstart'] > row['send']) + 1

	alignment = Alignment(
		name.replace('_', ' '), part,
		pd.Interval(row['sstart'] * modi, row['send'] * modi),
		int(row['qcovs']), row['bitscore']
	)

	sequences[row['sseqid']].alignments.append(alignment)

parts_order = {}
parts_required = set()
for o, part in enumerate(args.parts.split(',')):
	if '(o)' in part:
		part = part.replace('(o)', '').strip()
	else:
		parts_required.add(part)
	
	parts_order[part.strip()] = o

constructs = []
parts = []
for sequence in sequences.values():
	filtered = []
	parts_filtered = []

	sequence.alignments.sort(key=lambda x: (x.score, x.coverage), reverse=True)
	for alignment in sequence.alignments:
		olaps = np.array([overlap(alignment.interval, a.interval) for a in filtered])
		if not any(olaps > (0.5 * alignment.interval.length)):
			parts_filtered += alignment.part.split('-')
			filtered.append(alignment)

			sequence.seq.features.append(SeqFeature(makeFeatureLocation(alignment), id=alignment.name, type="gene", qualifiers={'label': alignment.name}))

	missing_parts = len(parts_required - set(parts_filtered))

	filtered.sort(key=lambda a: a.interval.left)
	for alignment in filtered:
		parts.append((sequence.name, missing_parts, alignment.coverage, alignment.score, alignment.part, alignment.name))

	filtered = reorderAlignments(filtered, parts_order)
	constructs.append((sequence.name, len(sequence.seq), missing_parts,
		'|'.join([a.name for a in filtered]))
	)

if args.annotate:
	with ZipFile(args.output / f'annotated.zip', 'w') as zip:
		for sequence in sequences.values():
			with zip.open(f'{sequence.name}.gb', 'w') as f:
				SeqIO.write(sequence.seq, TextIOWrapper(f), "gb")

constructs = pd.DataFrame(constructs, columns=('read', 'length', 'missing_parts', 'name'))
constructs.to_csv(args.output / f'constructs.csv')

num_constructs = len(pd.unique(constructs[constructs.name != ''].read))
complete = len(pd.unique(constructs[constructs.missing_parts < 1].read))

parts = pd.DataFrame(parts, columns=('read', 'missing_parts', 'coverage', 'score', 'type', 'name'))
parts.to_csv(args.output / f'parts.csv')

print(f'Found {num_constructs} constructs ({complete} complete) in {len(sequences)} sequences.')

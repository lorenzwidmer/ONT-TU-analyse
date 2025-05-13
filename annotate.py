import subprocess
import argparse
import tempfile
import shutil
import itertools

from pathlib import Path
from collections import namedtuple
from zipfile import ZipFile
from io import TextIOWrapper

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

import pandas as pd
import numpy as np

Sequence = namedtuple('Sequence', ['name', 'seq', 'alignments'])
Alignment = namedtuple('Alignment', ['name', 'part', 'positions', 'posnums', 'interval', 'direction', 'coverage', 'score'])

parser = argparse.ArgumentParser(description='''Script to find muCAR receptors in NanoPore Sequencing files.''', 
	prog='python3 -m sbo-builder')
parser.add_argument('fastqs', nargs='+', metavar='seq.fastq', type=Path,
	help='File containing the raw sequencing data.')
parser.add_argument('-o', '--output', type=Path, default=Path('./'),
	help='Name of the output file.')
parser.add_argument('-a', '--annotate', dest='annotate', action='store_true',
	help='Create zip file with all annotated sequences.')
parser.add_argument('-b', '--bidir', dest='annotate', action='store_true',
	help='Bi-directional reads.')
parser.add_argument('-p', '--parts', dest='parts', type=str, default='promo,rbs,signal,ntag,cds,ctag,stop,term',
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

def splitTUs(alignment, jump=-1):
	positions = [a.positions for a in alignment]
	positions = np.array(list(itertools.chain.from_iterable(positions)))

	splits = (np.argwhere(np.diff(positions) < jump).flatten() + 1).tolist()
	splits = [0] + splits + [len(alignment)]

	return [alignment[splits[i - 1]:splits[i]] for i in range(1, len(splits))]

root = Path(__file__).parent
folder = root / 'blast'
args.output.mkdir(exist_ok=True)

headers = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qcovs']
part_order = {p.strip():i for i, p in enumerate(args.parts.split(','))}

sequences = {}
for fastq in args.fastqs:
	with fastq.open('r') as f:
		fasta_sequences = SeqIO.parse(f, 'fastq')

		for fasta in fasta_sequences:
			fasta.id = f'{fasta.id}_{len(sequences)}'
			fasta.annotations["molecule_type"] = "DNA"
			sequences[fasta.id] = Sequence(fasta.id, fasta, [])

if not (args.output / 'alignments.csv').is_file():
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
else:
	results = pd.read_csv(args.output / 'alignments.csv', sep="\t", names=headers)

for i, row in results.iterrows():
	name, part = row['qseqid'].split('|')
	positions = [part_order.get(x.strip(), 0) for x in part.split('-')]

	if row['sstart'] > row['send']:
		interval = pd.Interval(row['send'], row['sstart'])
	else:
		interval = pd.Interval(row['sstart'], row['send'])

	alignment = Alignment(
		name.replace('_', ' '), part, positions,
		''.join([str(p) for p in positions]), interval,
		row['sstart'] > row['send'],
		int(row['qcovs']), row['bitscore']
	)

	sequences[row['sseqid']].alignments.append(alignment)

constructs = []
parts = []
for sequence in sequences.values():
	filtered = []

	sequence.alignments.sort(key=lambda x: (x.score, x.coverage), reverse=True)
	for alignment in sequence.alignments:
		olaps = np.array([overlap(alignment.interval, a.interval) for a in filtered])
		if not any(olaps > (0.5 * alignment.interval.length)):
			filtered.append(alignment)

			sequence.seq.features.append(SeqFeature(makeFeatureLocation(alignment), id=alignment.name, type="gene", qualifiers={'label': alignment.name}))

	filtered.sort(key=lambda a: a.interval.left)

	splitty = splitTUs(filtered)
	for tu, split in enumerate(splitty):
		for alignment in split:
			parts.append((sequence.name, f'TU{tu+1}', alignment.coverage, alignment.score, alignment.part, alignment.name))

		constructs.append((sequence.name, f'TU{tu+1}', len(sequence.seq),
			''.join([f.posnums for f in split]),
			'|'.join([f.name for f in split])
		))

if args.annotate:
	with ZipFile(args.output / f'annotated.zip', 'w') as zip:
		for sequence in sequences.values():
			with zip.open(f'{sequence.name}.gb', 'w') as f:
				SeqIO.write(sequence.seq, TextIOWrapper(f), "gb")

constructs = pd.DataFrame(constructs, columns=('read', 'TU', 'length', 'parts', 'name'))
constructs.to_csv(args.output / f'constructs.csv')

parts = pd.DataFrame(parts, columns=('read', 'TU', 'coverage', 'score', 'type', 'name'))
parts.to_csv(args.output / f'parts.csv')

print(f'Found {len(constructs)} constructs in {len(sequences)} sequences.')

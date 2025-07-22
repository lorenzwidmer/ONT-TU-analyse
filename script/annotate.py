import subprocess
import argparse
import tempfile
import shutil
import json

from pathlib import Path
from zipfile import ZipFile
from io import TextIOWrapper

from Bio import SeqIO

import pandas as pd

from classes import Sequence, Alignment 

parser = argparse.ArgumentParser(description='''Script to annotate Microsynth full plasmid seq ONT reads. The minimal input is a parts.fasta file with a list of possible sequences in the library and a clean raw-reads.fastq file.''', 
	prog='python3 annotate.py')
parser.add_argument('partfile', metavar='PART', type=Path,
	help='File containing the DNA parts (.fasta).')
parser.add_argument('fastqs', nargs='+', metavar='SEQ', type=Path,
	help='File containing the raw sequencing data.')
parser.add_argument('-o', '--output', type=Path, default=Path('./'),
	help='Name of the output file.')
parser.add_argument('-a', '--annotate', dest='annotate', action='store_true',
	help='Create zip file with all annotated sequences.')
parser.add_argument('-b', '--bidir', dest='annotate', action='store_true',
	help='Bi-directional reads.')
parser.add_argument('-p', '--parts', dest='parts', type=str, default='promo,rbs,signal,ntag,cds,ctag,stop,term',
	help='Part categories and their order in the TU.')
parser.add_argument('-f', '--filter', type=str, default='{"TU1":[0,1,2,5], "TU2":[0,1,2,4]}',
	help='Part filter applied to the parts in the cnstruct.')

args = parser.parse_args()

root = Path(__file__).parent
folder = root / 'blast'
args.output.mkdir(exist_ok=True)

headers = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qcovs']
part_order = {p.strip():i for i, p in enumerate(args.parts.split(','))}

part_filter = {k: set(f) for k, f in json.loads(args.filter).items()}

sequences = {}
for fastq in args.fastqs:
	with fastq.open('r') as f:
		fasta_sequences = SeqIO.parse(f, 'fastq')

		for fasta in fasta_sequences:
			fasta.id = f'{fasta.id}_{len(fasta.seq)}'
			fasta.annotations["molecule_type"] = "DNA"
			sequences[fasta.id] = Sequence(fasta.id, fasta)

if not (args.output / 'alignments.csv').is_file():
	with tempfile.TemporaryDirectory() as tempdir:
		tempdir = Path(tempdir)

		with (tempdir / 'plasmids.fasta').open("w") as f:
			SeqIO.write([sequence.sequence for sequence in sequences.values()], f, "fasta")

		if not (tempdir / 'plasmids.ndb').is_file():
			subprocess.run([
				'makeblastdb', '-in', 'plasmids.fasta', '-dbtype', 'nucl', '-parse_seqids', '-out', 'plasmids', '-title', 'no-idea'
			], cwd=tempdir)

			subprocess.run([
				'blastn', '-out', tempdir / 'alignments.csv', '-outfmt', ' '.join(['6'] + headers),
				'-max_target_seqs', str(len(sequences)),
				'-query', args.partfile.resolve(), '-db', tempdir / 'plasmids'
			], cwd=tempdir)

		results = pd.read_csv(tempdir / 'alignments.csv', sep="\t", names=headers)
		shutil.copy(tempdir / 'alignments.csv', args.output / 'alignments.csv')
else:
	results = pd.read_csv(args.output / 'alignments.csv', sep="\t", names=headers)

for i, row in results.iterrows():
	sequences[row['sseqid']].add(Alignment(row, part_order))

constructs = []
parts = []
for sequence in sequences.values():
	sequence.filter()

	splitty = sequence.splitTUs()
	for tu, split in enumerate(splitty):
		for alignment in split.alignments:
			for part in alignment.parts:
				parts.append((split.name, f'TU{tu+1}', alignment.coverage, alignment.score, part.name, alignment.name))

		construct = split.partDict(filter=part_filter.get(f'TU{tu+1}', None))
		constructs.append((split.name, f'TU{tu+1}', len(split.sequence),
			''.join([str(n) for n, _, _ in construct]),
			'|'.join([n for _, _, n in construct])
		))

if args.annotate:
	with ZipFile(args.output / f'annotated.zip', 'w') as zip:
		for sequence in sequences.values():
			with zip.open(f'{sequence.name}.gb', 'w') as f:
				SeqIO.write(sequence.sequence, TextIOWrapper(f), "gb")

constructs = pd.DataFrame(constructs, columns=('read', 'TU', 'length', 'parts', 'name'))
constructs.to_csv(args.output / f'constructs.csv')

parts = pd.DataFrame(parts, columns=('read', 'TU', 'coverage', 'score', 'type', 'name'))
parts.to_csv(args.output / f'parts.csv')

print(f'Found {len(constructs)} constructs in {len(sequences)} sequences.')

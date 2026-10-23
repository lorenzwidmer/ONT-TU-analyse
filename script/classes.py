import itertools
import gzip
import subprocess
import multiprocessing
import shutil

from pathlib import Path

from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

import pandas as pd
import numpy as np

class BlastSearch:
    HEADERS = ['qseqid', 'sseqid', 'pident', 'length', 'slen', 'mismatch', 'gapopen', 'qstart', 'qend', 'evalue', 'bitscore'] 

    def __init__(self, title, folder):
        self.title = title
        self.folder = Path(folder)

        self.seqfile = self.folder / f'blast-{self.title}.fasta'
        self.dbfile = self.folder / f'blast-parts-{self.title}.ndb'
        self.resultfile = self.folder / f'blast-raw-{self.title}.csv'

    def createDB(self, sequences):
        with self.seqfile.open("w") as f:
            SeqIO.write(sequences, f, "fasta")

        subprocess.run([
            'makeblastdb', '-in', self.seqfile.resolve(), '-dbtype', 'nucl', '-parse_seqids', '-out', self.dbfile.stem, '-title', self.title
        ], cwd=self.folder)
    
    def runBLAST(self, queryfile, numSeqs, task='blastn'):
        subprocess.run([
            'blastn', '-out', self.folder / f'blast-raw-{self.title}.csv', '-outfmt', ' '.join(['6'] + BlastSearch.HEADERS),
            '-max_target_seqs', str(numSeqs),
            '-num_threads', str(multiprocessing.cpu_count()),
            '-task', task,
            '-query', queryfile.resolve(), '-db', self.folder / self.dbfile.stem
        ], cwd=self.folder)

    def run(self, sequences, queryfile, task='blastn', csvout=None, force=False):
        if not self.dbfile.is_file():
            self.createDB(sequences)
        if force or not self.resultfile.is_file():
            self.runBLAST(queryfile, len(sequences), task)
        if csvout is not None:
            shutil.copy(self.resultfile, csvout)

        return pd.read_csv(self.resultfile, sep="\t", names=BlastSearch.HEADERS)

class FastxFile:
    def __init__(self, path):
        self.path = Path(path)

    def __enter__(self):
        with self.path.open('rb') as f:
            self.gzip = f.read(2) == b'\x1f\x8b'

        if self.gzip:
            self.f = gzip.open(self.path, "rt")
        else:
            self.f = self.path.open('r')
        return self.f

    def __exit__(self, type, value, traceback):
        self.f.close()


class Part:
    def __init__(self, name, part_order={}):
        self.name = name
        self.number = part_order.get(name.strip(), 0)

class Alignment:
    def __init__(self, row, part_order={}):
        name, part = row['sseqid'].replace('gnl|', '').split('|')
        
        self.name = name.replace('_', ' ')
        self.part = part
        self.parts = [Part(p, part_order) for p in part.split('-')]
        
        self.direction = row['qstart'] > row['qend']
        if row['qstart'] > row['qend']:
            self.interval = pd.Interval(row['qend'], row['qstart'])
        else:
            self.interval = pd.Interval(row['qstart'], row['qend'])
        
        self.coverage = int(row['scovs'])
        self.score = row['bitscore']

    def merge(self, other):
        self.interval = pd.Interval(min(self.interval.left, other.interval.left), max(self.interval.right, other.interval.right))
        self.coverage = self.coverage + other.coverage
        self.score = (self.score + other.score) / 2

def overlap(a, b):
	return min(a.right, b.right) - max(a.left, b.left) + 1

class Sequence:
    def __init__(self, name, sequence):
        self.alignments = []
        self.sequence = sequence
        self.name = name

    @property
    def parts(self):
        return list(itertools.chain.from_iterable([a.parts for a in self.alignments]))

    def add(self, alignment):
        self.alignments.append(alignment)

    def filter(self, merge_same=True):
        filtered = []
        alignments = sorted(self.alignments, key=lambda a: (a.score, a.coverage), reverse=True)

        for a in alignments:
            olaps = np.array([overlap(a.interval, f.interval) for f in filtered])
            if not any(olaps > (0.5 * a.interval.length)):
                filtered.append(a)

                self.sequence.features.append(
                    SeqFeature(FeatureLocation(a.interval.left, a.interval.right, strand=a.direction),
                        id=a.name, type="gene",
                        qualifiers={'label': a.name}
                ))
        
        self.alignments = sorted(filtered, key=lambda a: a.interval.left)

        if merge_same and len(self.alignments) > 0:
            filtered = [self.alignments[0]]
            for i in range(1, len(self.alignments)):
                firsty = self.alignments[i-1]
                secondy = self.alignments[i]

                if firsty.name == secondy.name and firsty.part == secondy.part and firsty.direction == secondy.direction:
                    firsty.merge(secondy)
                else:
                    filtered.append(secondy)
            self.alignments = filtered

    def splitTUs(self, jump=-1):
        positions = np.array([(a.parts[0].number, a.parts[-1].number) for a in self.alignments]).flatten()

        splits = ((np.argwhere(np.diff(positions) <= jump).flatten() + 1) / 2).tolist()
        splits = [0] + splits + [len(self.alignments)]

        for i in range(1, len(splits)):
            seq = Sequence(self.name, self.sequence)
            seq.alignments = self.alignments[int(splits[i - 1]):int(splits[i])]
            yield seq

    def partDict(self, filter=None):
        output = []
        for alignment in self.alignments:
            for p in alignment.parts:
                if filter is not None and p.number in filter:
                    output.append((p.number, p.name, alignment.name))

        return output


                

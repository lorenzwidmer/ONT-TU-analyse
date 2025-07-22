import itertools

from collections import namedtuple
from Bio.SeqFeature import SeqFeature, FeatureLocation

import pandas as pd
import numpy as np

Part = namedtuple('Part', ['name', 'number'])

class Part:
    def __init__(self, name, part_order={}):
        self.name = name
        self.number = part_order.get(name.strip(), 0)

class Alignment:
    def __init__(self, row, part_order={}):
        name, part = row['qseqid'].split('|')
        
        self.name = name.replace('_', ' ')
        self.part = part
        self.parts = [Part(p, part_order) for p in part.split('-')]
        
        self.direction = row['sstart'] > row['send']
        if row['sstart'] > row['send']:
            self.interval = pd.Interval(row['send'], row['sstart'])
        else:
            self.interval = pd.Interval(row['sstart'], row['send'])
        
        self.coverage = int(row['qcovs'])
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
        positions = np.array([a.parts[0].number for a in self.alignments])

        splits = (np.argwhere(np.diff(positions) <= jump).flatten() + 1).tolist()
        splits = [0] + splits + [len(self.alignments)]

        for i in range(1, len(splits)):
            seq = Sequence(self.name, self.sequence)
            seq.alignments = self.alignments[splits[i - 1]:splits[i]]
            yield seq

    def partDict(self, filter=None):
        output = []
        for alignment in self.alignments:
            for p in alignment.parts:
                if filter is not None and p.number in filter:
                    output.append((p.number, p.name, alignment.name))

        return output


                

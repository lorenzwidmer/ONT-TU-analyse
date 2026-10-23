import csv
import re
import argparse

from typing import Any, Dict

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

@dataclass
class Variant:
	name: str
	gaps: int
	count: int
	count2: float

	def reverse(self, structure):
		return [part.getName(self.name[part.index]) for part in structure.values()]
	
@dataclass
class Part:
	tu: str
	type: str
	index: int
	bindings: Dict[str, str]
	revbindings: Dict[str, str]
	combinable: bool

	def getBinding(self, name):
		if name not in self.bindings:
			binding = chr(65 + len(self.bindings))
			self.bindings[name] = binding
			self.revbindings[binding] = name
		return self.bindings[name]
	
	def getName(self, binding):
		return self.revbindings.get(binding, None)

	
def transformStructure(structure):
	output = {}
	for tu in structure.split(';'):
		tu, parts = tu.split(':')
		tu = tu.strip()
		
		for part in parts.split(','):
			combinable = not part.strip().startswith('!')
			part = part.strip(' !')
			output[f'{tu}.{part}'] = Part(tu, part, len(output), {}, {}, combinable)
	
	return output

def count(constructs, structure):
	output = {}
	mask = [part.index for part in structure.values() if not part.combinable]

	for c in constructs.values():
		for x in mask:
			if c[x] == '.':
				c[x] = 'X'

		gaps = c.count('.')

		construct = ''.join(c)
		if construct not in output:
			output[construct] = Variant(construct, gaps, 1, 0)
		else:
			output[construct].count += 1
	
	return output

parser = argparse.ArgumentParser(prog='PROG', usage='%(prog)s [options]')
parser.add_argument('files', nargs='+', type=Path, help='bar help')
parser.add_argument('-s', '--structure', default='TU1: promo,rbs,signal,cds,ctag; TU2: promo,rbs,signal,!cds,ctag', type=str, help='bar help')
parser.add_argument('-o', '--out', default=Path('./out.csv'), type=Path, help='bar help')
parser.add_argument('-g', '--gaps', default=2, type=int, help='bar help')
args = parser.parse_args()

structure = transformStructure(args.structure)

stats = []
output = []

for file in args.files:
	constructs = defaultdict(lambda: ['.']*len(structure))
	with file.open('r', newline='') as csvfile:
		spamreader = csv.reader(csvfile, delimiter=',', quotechar='"')
		for row in spamreader:
			_, read, tu, _, _, typ, name = row
			idx = f'{tu}.{typ}'

			if idx in structure:
				part = structure[idx]
				constructs[read][part.index] = part.getBinding(name)

	
	counts = count(constructs, structure)
	fullcounts = {k: c for k, c in counts.items() if c.gaps == 0}
	fullconstructs = list(fullcounts.keys())
	for _, construct in counts.items():
		if construct.gaps > 0 and construct.gaps <= args.gaps:
			needle = re.compile(f'^{construct.name}$')
			matches = [s for s in fullconstructs if needle.match(s)]

			if len(matches) > 0:
				totally = sum([fullcounts[match].count for match in matches])
				for match in matches:
					fullcounts[match].count2 += (fullcounts[match].count/totally) * construct.count
			elif construct.count > 1:
				fullcounts[construct.name] = construct
				
		elif construct.gaps > 0 and construct.count > 1:
			fullcounts[construct.name] = construct

	stat = pd.DataFrame([(c.gaps, c.count) for c in counts.values()], columns=["gaps", "count"])
	stat = stat.groupby(['gaps'], sort=True).sum()
	stats.append([file.name] + list(stat['count']))

	for c in sorted(fullcounts.values(), key=lambda c: c.count + c.count2, reverse=True):
		output.append([
			file.name, c.name, c.gaps, c.count, c.count2
		] + c.reverse(structure))

print(stats)
output = pd.DataFrame(output, columns=['file', 'name', 'gaps', 'count1', 'count2'] + list(structure.keys()))
output.to_csv(args.out)



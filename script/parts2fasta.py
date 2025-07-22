
import re
import csv

import argparse

from pathlib import Path

parser = argparse.ArgumentParser(description='''Script to convert a parts.csv benchling export to a fasta file that can be read by annotate.py.''', 
	prog='python3 parts2fasta.py')
parser.add_argument('input', metavar='IN', type=Path,
	help='Benchling export (.csv).')
parser.add_argument('output', metavar='OUT', type=Path,
	help='Output file (.fasta).')

args = parser.parse_args()

bsaI_regex = re.compile('GGTCTC[ACGT]([ACGT]{4})([ACGT]+)([ACGT]{4})[ACGT]GAGACC', flags=re.IGNORECASE)
lvl0 = ["CTAT", "ACTT", "CATA", "GGAA", "AGTG", "ACCG", "GGCT", "CGAC", "TGTT"]
lvl0_dict = {overhang:i for i,overhang in enumerate(lvl0)}
partnames = ['promo', 'rbs', 'signal', 'ntag', 'cds', 'ctag', 'stop', 'term']
output = []

with args.input.open(newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',', quotechar='"')
    for row in spamreader:
        fullname, author, seq, lvl = row[0:4]

        if lvl.startswith('lvl') and lvl != 'lvl0':
            continue

        match = bsaI_regex.search(seq.upper())
        if match is not None:
            start = lvl0_dict.get(match.group(1), 0)
            end = lvl0_dict.get(match.group(3), 0)

            name = fullname.strip('✅☑️ ').replace(' ', '-').replace(' ', '-')
            name = bytes(name, 'utf-8').decode('utf-8', 'ignore')

            parts = '-'.join(partnames[start:end])

            output.append(f">{name}|{parts}\n")
            output.append(f"{match.group(2)}\n")
        else:
            print(f'{fullname} contains no bsai')

with args.output.open('w') as f:
    f.writelines(output)

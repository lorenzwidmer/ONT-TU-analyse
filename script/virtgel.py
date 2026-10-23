import sys
import argparse
from pathlib import Path

import itertools
import pysam
import numpy as np
import plotly.express as px

parser = argparse.ArgumentParser(prog='PROG', usage='%(prog)s [options]')
parser.add_argument('files', nargs='+', type=Path, help='bar help')
parser.add_argument('-r', '--resolution', type=int, help='Number of bins', default=100)
parser.add_argument('-l', '--min_length', type=int, help='Lower length cutoff', default=0)
parser.add_argument('-L', '--max_length', type=int, help='Upper length cutoff', default=None)
parser.add_argument('-n', '--no_normalise', action='store_true', help='Upper length cutoff')
args = parser.parse_args()

# Collect read lengths
lengths = {}

for file in args.files:
	if file.suffix == '.bam':
		with pysam.AlignmentFile(file, "rb") as bam:
			lengths[file.stem] = [read.query_length for read in bam if not read.is_unmapped]
	else:
		with pysam.FastxFile(file) as fastx:
			lengths[file.stem] = [len(read.sequence) for read in fastx]

lenghtLim = [args.min_length, args.max_length]
if lenghtLim[1] is None:
	abc = np.array([l for l in itertools.chain(*lengths.values())])
	lenghtLim[1] = np.quantile(abc, 0.95) * 1.1

maxlen = max(map(max, lengths.values()))
bins = np.arange(lenghtLim[0], lenghtLim [1], (lenghtLim[1] - lenghtLim[0]) / args.resolution)

heatmap = np.array([
    np.histogram(lengths, bins=bins)[0]
    for lengths in lengths.values()
])

if not args.no_normalise:
	heatmap = 100 * heatmap / np.sum(heatmap, axis=1, keepdims=True)

fig = px.imshow(
    heatmap.T,
    x=list(lengths.keys()),
    y=bins[:-1],
    labels={"x": "Sample", "y": "Length", "color": "Percentage"},
	color_continuous_scale='deep'
)

fig.write_html("plot.html")
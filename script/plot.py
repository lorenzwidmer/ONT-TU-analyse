import argparse
import re
import sys

from pathlib import Path

import pandas as pd
import numpy as np

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def parseGroupBy(groupby):
	output = {}
	for tu in groupby.split(';'):
		if ':' not in tu:
			continue

		tu, parts = tu.split(':')
		tu = tu.strip()
		
		output[tu] = [f'{tu}.{part.strip()}' for part in parts.split(',')]

	return output

parser = argparse.ArgumentParser(prog='PROG', usage='%(prog)s [options]')
parser.add_argument('samples', type=Path, help='bar help')
parser.add_argument('layouts', type=Path, help='bar help')
parser.add_argument('counts', type=Path, nargs='+', help='bar help')

parser.add_argument('--sorts', default=None, help='bar help')
parser.add_argument('--no-partial', dest="no_partial", action='store_true',
				help='Upper length cutoff')
args = parser.parse_args()

# Load CSV
df = pd.concat([pd.read_csv(c) for c in args.counts])
samples = pd.read_csv(args.samples)


df = pd.merge(df, samples, on="file")

if args.no_partial:
	df["countA"] = df['count1']
else:
	df["countA"] = df['count1'] + df['count2']

if args.sorts is not None:
	sorts = set([int(sort) for sort in args.sorts.split(',')])
	df = df.query("sort in @sorts")

layouts = pd.read_csv(args.layouts, comment='#')
layouts = layouts.loc[layouts['row'] >= 0]

layout_setting = {}
with Path(args.layouts).open('r') as f:
	layout_cmd = next(f, '#').strip("# ")
	layout_setting = {s[0].strip(): s[1].strip() for s in re.findall(r'([^;]+)=\s*"([^;]+)"\s*;', layout_cmd)}

facetcol = layout_setting.get('facetwrap', 'TU1.cds')
colorcol = layout_setting.get('color', facetcol)
min_scount = int(layout_setting.get('min_scount', 50))
colors = ['#FFAD0A', '#E85D75', '#25CED1', '#E086D3', '#462749']

facets = {facet: i+1 for i, facet in enumerate(np.sort(df[facetcol].dropna().unique()))}
colors = {color: colors[i] for i, color in enumerate(np.sort(df[colorcol].dropna().unique()))}

fig = make_subplots(
	rows=layouts.shape[0],
	cols=len(facets),
	vertical_spacing=0.02,
	horizontal_spacing=0.02,
	shared_xaxes=True,
	shared_yaxes=True,
	column_titles=list(facets.keys()) ,
	row_titles=list(layouts['title'])
)

for plot_row, layout in enumerate(layouts.sort_values(by=["row"]).iterrows()):
	_, layout = layout

	subdf = df.query(f'groupname == "{layout.loc['sample']}" and ({layout.loc['filter']})').copy()
	if subdf.shape[0] <= 0:
		continue

	groupby = parseGroupBy(layout.loc['groupby'])
	for tu, cols in groupby.items():
		subdf[tu] = subdf[cols].apply(lambda s: s.str.cat(sep=', '), axis=1)
	tus = list(groupby.keys())
	
	group = subdf.sort_values("sort")

	sums = group.groupby(['sort'])["countA"].transform("sum")
	group['countN'] = 100 * group["countA"] / sums

	#filter sorts with total counts < min_scount
	group = group[group.groupby("sort")["countA"].transform("sum") > min_scount]

	for name, subdfy in group.groupby(tus):
		at = subdfy['TU1.ctag'].iloc[0]
		signal = subdfy["TU1.signal"].iloc[0]
		if at not in facets:
			continue

		trace = subdfy.groupby("sort")["countN"].sum()
		fig.add_trace(
			go.Scatter(
				x=list(trace.index),
				y=list(trace),
				name=str(subdfy["name"].iloc[0]),
				mode="lines+markers",
				showlegend=False,
				customdata=subdfy[tus],
				line=dict(
        			color=colors.get(signal, "black"),
					width=1.5
				),
				marker=dict(size=4)
			),
			col=facets[at],
			row=plot_row + 1
		)

	sums = group.groupby(['sort']).agg({'countA': 'sum'})
	offset = max(20, group.groupby(tus + ['sort'])["countN"].sum().max())
	for idx, row in sums.iterrows():
		for coly in facets.values():
			fig.add_annotation(x=idx, y=-0.1 * offset,
				text=round(row.loc['countA']),
				showarrow=False,
				yshift=0,
				col=coly,
				row=plot_row + 1)

	# Show x-axis labels as S0, S1, S2, ...
	sort_values = sorted(group["sort"].unique())
	fig.update_xaxes(
		tickvals=sort_values,
		ticktext=[f"S{s}" for s in sort_values]
	)

	fig.update_yaxes(range=[-0.2 * offset, offset * 1.1], row=plot_row + 1)
	fig.update_layout(template="simple_white")

	abc = '<br>'.join([f"{tu}: <i>%{{customdata[{i}]}}</i>" for i, tu in enumerate(tus)])
	fig.update_traces(
		hovertemplate = f"""<b>%{{y:.2f}} %</b><br>{abc}
		"""
	)

with open('p_graph.html', 'w') as f:
	f.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))
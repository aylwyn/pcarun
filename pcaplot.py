#!/usr/bin/env python
# Aylwyn Scally 2014

import numpy as np
import matplotlib.pyplot as plt
#import glob
import pandas as pd
import argparse
#from scipy.stats import norm
from matplotlib import rcParams

class PCAGroup(object):
	def __init__(self, name, col='blue', ptype = 'o'):
		self.name = name
		self.pcol = col
		self.ptype = ptype

	def show(self):
		print("\t".join([self.name, self.ptype, self.pcol]))

p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#pp.add_argument('--sim', action='store_true', default = False, help = 'dry run')
p.add_argument('--hc', action='store_true', default = False, help = 'hard copy (pdf)')
p.add_argument('--labels', action='store_true', default = False, help = 'label points')
p.add_argument('--eigenstrat', action='store_true', default = False, help = 'eigenstrat format')
p.add_argument('--opensymbols', action='store_true', default = False, help = 'plot open symbols')
p.add_argument('--skip', type=int, default = 0, help = 'skip initial lines')
p.add_argument('EVEC_FILE')
p.add_argument('-c', '--components', default = '1,2', help = 'Components to plot')
p.add_argument('-L', '--legend_file', default = '', help = 'File with info for legend')
p.add_argument('--nolegend', action='store_true')
p.add_argument('-e', '--eigenval_file', default = '', help = 'File of eigenvalues')
p.add_argument('-o', '--out', default = '', help = 'prefix for pdf file')
p.add_argument('--figsize', default = '2.25,2.25', help = 'figure size in inches')
p.add_argument('--fontsize', default = 6, help = 'figure font size in points')
args = p.parse_args()

if args.hc:
	rcParams['axes.labelsize'] = 6
	rcParams['xtick.labelsize'] = 5
	rcParams['ytick.labelsize'] = 5
	rcParams['legend.fontsize'] = 5
fsize = [float(x) for x in args.figsize.split(',')]

msize = 5

fig, ax = plt.subplots(1, 1, squeeze=False)
axv = ax.reshape(-1)

plt.sca(axv[0])

if args.eigenstrat:
	namecol = 0
else:
	namecol = 1
pc = [int(x) for x in args.components.split(',')]

pgrp = {}
group = {}
if args.legend_file:
	for line in open(args.legend_file):
		if line.isspace():
			continue
		tok = line.strip().split('\t')
		if tok[0] == 'group':
			group[tok[1]] = PCAGroup(*tok[1:])
#			group[tok[1]].show()
		else:
			pgrp[tok[0]] = group[tok[1]]

tb = pd.read_table(args.EVEC_FILE, header = None, comment = '#', sep = '\s*', skiprows = args.skip)

usedgrps = set()
for i in range(len(tb)):
	name = tb.ix[i, namecol]
	if name in pgrp:
		tpinfo = pgrp[name]
	else:
		tpinfo = PCAGroup('')
#	print(pgrp[name].name)
	usedgrps.add(pgrp[name].name)
	if args.opensymbols:
		plt.plot(tb.ix[i, pc[0] + 1], tb.ix[i, pc[1] + 1], marker=tpinfo.ptype, markeredgecolor=tpinfo.pcol, markerfacecolor='none', markersize=msize)
	else:
		plt.plot(tb.ix[i, pc[0] + 1], tb.ix[i, pc[1] + 1], marker=tpinfo.ptype, color=tpinfo.pcol)
	if args.labels:
		plt.text(tb.ix[i, pc[0] + 1], tb.ix[i, pc[1] + 1], name, color='grey', fontsize='xx-small')

if args.eigenval_file:
	evals = [float(x) for x in open(args.eigenval_file).readlines()]
	evalsum = sum(evals)
	varfrac = [evals[x - 1] / evalsum for x in pc]
	plt.xlabel('PC%s (%.1f %% of variance)' % (str(pc[0]), 100 * varfrac[0]), labelpad=3)
	plt.ylabel('PC%s (%.1f %% of variance)' % (str(pc[1]), 100 * varfrac[1]), labelpad=2)
else:
	plt.xlabel('PC' + str(pc[0]), labelpad=3)
	plt.ylabel('PC' + str(pc[1]), labelpad=2)

#plt.xlim(0.05,0.15)
plt.axis([1.1*x for x in plt.axis()])

##plt.subplots_adjust(wspace = 0.1)
#fig = plt.gcf()
#fig.set_size_inches(12, 4)

lugrps = sorted(list(usedgrps))
if args.legend_file and not args.nolegend:
	from matplotlib.lines import Line2D # Using Line2D as proxy artist
	arts = []
#	for g in group.keys():
	for g in lugrps:
		if args.opensymbols:
			arts.append(Line2D(range(1), range(1), color="white", marker=group[g].ptype, markeredgecolor=group[g].pcol, markerfacecolor='none', markersize=msize))
		else:
			arts.append(Line2D(range(1), range(1), color="white", marker=group[g].ptype, markerfacecolor=group[g].pcol))
	 
	plt.legend(arts, lugrps, loc = "best", numpoints=1)#, fontsize='12')

plt.subplots_adjust(left=0.36/fsize[0],bottom=0.26/fsize[1], right=1 - 0.04/fsize[0], top= 1 - 0.04/fsize[1], hspace = 0.25)

if args.hc:
	fig.set_size_inches(fsize[0], fsize[1])
	if not args.out:
		args.out = args.EVEC_FILE
	ofile = args.out + '.PC' + '-'.join([str(x) for x in pc]) + '.pdf'
	plt.savefig(ofile)
	print('plot written to %s' % ofile)
else:
	plt.show()

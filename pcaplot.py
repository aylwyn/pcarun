#!/usr/bin/env python
# Aylwyn Scally 2014

import numpy as np
import matplotlib.pyplot as plt
#import glob
import pandas as pd
import argparse
#from scipy.stats import norm

class PCAGroup(object):
    def __init__(self, name, col='blue', ptype = 'o'):
        self.name = name
        self.pcol = col
        self.ptype = ptype

    def show(self):
        print("\t".join([self.name, self.ptype, self.pcol]))

p = argparse.ArgumentParser()
#pp.add_argument('--sim', action='store_true', default = False, help = 'dry run')
p.add_argument('--hc', action='store_true', default = False, help = 'hard copy (pdf)')
p.add_argument('--labels', action='store_true', default = False, help = 'label points')
p.add_argument('--eigenstrat', action='store_true', default = False, help = 'eigenstrat format')
p.add_argument('--skip', type=int, default = 0, help = 'skip initial lines')
p.add_argument('EVEC_FILE')
p.add_argument('-c', '--components', default = '1,2', help = 'Components to plot [\'1,2\']')
p.add_argument('-L', '--legend_file', default = '', help = 'File with info for legend')
p.add_argument('-o', '--out', default = '', help = 'prefix for pdf file')
args = p.parse_args()

if args.eigenstrat:
	namecol = 0
	pc = [int(x) for x in args.components.split(',')]
else:
	namecol = 1
	pc = [int(x) + 1 for x in args.components.split(',')]

pgrp = {}
group = {}
if args.legend_file:
	for tok in (line.split() for line in open(args.legend_file)):
		if tok[0] == 'group':
			group[tok[1]] = PCAGroup(*tok[1:])
			group[tok[1]].show()
		else:
			pgrp[tok[0]] = group[tok[1]]

tb = pd.read_table(args.EVEC_FILE, header = None, comment = '#', sep = '\s*', skiprows = args.skip)

for i in range(len(tb)):
	name = tb.ix[i, namecol]
	if name in pgrp:
		tpinfo = pgrp[name]
	else:
		tpinfo = PCAGroup['0']
	plt.plot(tb.ix[i, pc[0]], tb.ix[i, pc[1]], marker=tpinfo.ptype, color=tpinfo.pcol)
	if args.labels:
		plt.text(tb.ix[i, pc[0]], tb.ix[i, pc[1]], name, color='grey', fontsize='xx-small')

plt.xlabel('PC' + str(pc[0]))
plt.ylabel('PC' + str(pc[1]))

#nobs = len(oddsrat)
#meanodds = np.mean(oddsrat[:, 1])
#sdodds = np.std(oddsrat[:, 1])
#plt.hlines(meanodds, 0, xmax, color='m')
#plt.hlines(0, 0, xmax, color='k')
#siglevels = [norm.ppf(1-0.25/nobs, meanodds, sdodds), norm.ppf(0.25/nobs, meanodds, sdodds)]
#plt.hlines(siglevels, 0, xmax, color='m', linestyles='dashed')

#plt.vlines(gchrs['gpos'], -2, 2)
#plt.gca().set_xticks(gchrs['gmid'][:-1])
#plt.gca().set_xticklabels(gchrs['nname'][:-1], fontsize='small')

#plt.xlim(0.05,0.15)
plt.axis([1.1*x for x in plt.axis()])

##plt.subplots_adjust(wspace = 0.1)
#fig = plt.gcf()
#fig.set_size_inches(12, 4)

if args.legend_file:
	from matplotlib.lines import Line2D # Using Line2D as proxy artist
	arts = []
	for g in group.keys():
		arts.append(Line2D(range(1), range(1), color="white", marker=group[g].ptype, markerfacecolor=group[g].pcol))
	 
	plt.legend(arts, group.keys(), loc = "best", numpoints=1)

if args.hc:
	if not args.out:
		args.out = args.EVEC_FILE
	ofile = args.out + '.' + ''.join(pc) + '.pdf'
	plt.savefig(ofile)
	print('plot written to %s' % ofile)
else:
	plt.show()

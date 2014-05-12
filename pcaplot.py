#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
#import glob
import pandas as pd
import argparse
#from scipy.stats import norm

p = argparse.ArgumentParser()
#pp.add_argument('--sim', action='store_true', default = False, help = 'dry run')
p.add_argument('--hc', action='store_true', default = False, help = 'hard copy (pdf)')
p.add_argument('--labels', action='store_true', default = False, help = 'label points')
p.add_argument('--eigenstrat', action='store_true', default = False, help = 'eigenstrat format')
p.add_argument('--skip', type=int, default = 0, help = 'skip initial lines')
p.add_argument('EVEC_FILE')
p.add_argument('-c', '--components', default = '1,2', help = 'Components to plot [\'1,2\']')
p.add_argument('-l', '--legend_file', default = '', help = 'File with info for legend')
p.add_argument('-o', '--out', default = '', help = 'prefix for pdf file')
args = p.parse_args()

pc = args.components.split(',')

#coldict = {'gorberber': 'red', 'gorbergra': 'green', 'gorgorgor': 'blue', 'Gbb': 'red', 'Gbg': 'green', 'Ggg': 'blue', 'Ggd': 'magenta'}
#coldict = {'Gbb': 'red', 'Gbg': 'green', 'Ggg': 'blue', 'Ggd': 'magenta'}
#colnames = ['Name', 'PC1', 'PC2', 'PC3', 'PC4', 'species']
tb = pd.read_table(args.EVEC_FILE, header = None, names=colnames, comment = '#', sep = '\s*', skiprows = args.skip)

for i in range(len(tb)):
	sp = tb.ix[i, 'species']
	icol = coldict[sp]
	plt.plot(tb.ix[i, pc[0]], tb.ix[i, pc[1]], 'o', color=icol)
	if args.labels:
		plt.text(tb.ix[i, pc[0]], tb.ix[i, pc[1]], tb.ix[i, 'Name'], color='grey', fontsize='xx-small')
#	plt.plot(tb['PC1'], tb['PC2'], '.')

#plt.title(iname)
plt.xlabel(pc[0])
plt.ylabel(pc[1])

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

#plt.legend(tuple(coldict.keys()), tuple(coldict.values()))

from matplotlib.lines import Line2D # Using Line2D as proxy artist
arts = []
for sp in coldict.keys():
	arts.append(Line2D(range(1), range(1), color="white", marker='o', markerfacecolor=coldict[sp]))
 
 # Calling the handles and labels to create the legend, where the handles are the club and circle created previously, and the labels are what the markers are labeled in the legend. Also moves the legend outside the figure
plt.legend(arts, coldict.keys(), loc = "best", numpoints=1)#, bbox_to_anchor = (1, 0.5), numpoints = 1)

if args.hc:
	if not args.out:
		args.out = args.EVEC_FILE
	ofile = args.out + '.' + ''.join(pc) + '.pdf'
	plt.savefig(ofile)
	print('plot written to %s' % ofile)
else:
	plt.show()

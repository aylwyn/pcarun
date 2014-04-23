#!/software/bin/python
# Aylwyn Scally 2014

import sys
import getopt
#import os
import os.path
#from math import sqrt
import logging
from logging import error, warning, info, debug, critical

loglevel = logging.WARNING

#os.umask(0002)

def lout(fout, *args):
	fout.write('\t'.join([str(x) for x in args]) + '\n')

def usage():
	sys.stderr.write('usage: %s sample_info_file [varsites_file] [-p out_prefix]\n' % (os.path.basename(sys.argv[0])))

# defaults
header = False
outpref = 'eig'
oldsampfmt = False
appendout = False

try:
	opts, args = getopt.gnu_getopt(sys.argv[1:], 'p:', ['oldsampfmt'])
except getopt.GetoptError:
	usage()
	sys.exit(2)
for (oflag, oarg) in opts:
	if oflag == '-p':
		outpref = oarg
	if oflag == '--append':
		appendout = True
	if oflag == '--oldsampfmt':
		oldsampfmt = True

#if args:
#	sampfile = open(args[0])
#else:
#	usage()
if len(args) > 1:
	fin = open(args[1])
else:
	fin = sys.stdin
fout = sys.stdout

vcfsamps = fin.readline().split()[1:]
#indpos = dict([(x, i) for i, line in enumerate(varsamps)])
#print(indpos)

indfile = open(outpref + '.ind', 'w')
for samp in vcfsamps:
	if oldsampfmt:
		ind = samp.split('/')[-1].split('.')[0] # samples in vcf are mapdir/species_sex_name.bam
		indivdat = list(reversed(ind.split('_')))
	else: # samples in vcf are species-sex-name
		indivdat = list(reversed(samp.split('-')))
	indfile.write('\t'.join(indivdat) + '\n') # need to write name sex species for each indiv
#	else:
#		print('sample %s not in %s; excluding' % (ind, args[0]))

if appendout:
	writemode = 'a'
else:
	writemode = 'w'
snpfile = open(outpref + '.snp', writemode)
genofile = open(outpref + '.geno', writemode)

genodict = {'00':'2', '01':'1', '11':'0', '..':'9'}
for line in fin:
	tok = line.split()
	if tok[0] == 'X':
		tok[0] = '23'
	snpfile.write('\t'.join(['_'.join(['snp', tok[0], tok[1]]), tok[0], '0.0', tok[1]]) + '\n')
	genos = [genodict[x] for x in tok[2].split('-')]
	genofile.write(''.join(genos) + '\n')

snpfile.close()
genofile.close()

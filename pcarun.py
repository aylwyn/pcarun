#!/software/bin/python
# Aylwyn Scally 2014

import sys
import argparse
import os
import re
import os.path
import glob
import logging
from logging import error, warning, info, debug, critical

import aosutils

DATASUBDIR = 'data'

def prep(args):
	debug(args)
	# make dir
	if args.s1name == 'vcf':
		vcf_input = [(os.path.abspath(x), '') for x in args.VCF_FILE]
		if args.vcf_list:
			for line in open(args.vcf_list):
				vcf_input.append(tuple(line.split()))

		sample_input = []
		if args.samples:
			sample_input = args.samples.split(',')
		if args.sample_list:
			for line in open(args.sample_list):
				sample_input.append(line.strip())
		if not sample_input:
			sample_input = ['*']
		samp = ','.join(sample_input)

	elif args.s1name == 'ms':
		ms_file = os.path.abspath(args.MS_FILE)
		if not os.path.exists(ms_file):
			error('No such file: %s' % ms_file)
			return(2)

	debug(samp)
	if args.dirpref:
		sname = args.dirpref
	elif args.s1name == 'vcf':
		sname = os.path.basename(vcf_input[0][0]).replace('.bgz', '').replace('.vcf', '')
		warning('taking dirname from first vcf file')
	elif args.s1name == 'ms':
		sname = os.path.basename(ms_file)

	maindir = sname + '.pcadir'
	if not os.path.exists(maindir):
		info('creating directory %s' % maindir)
		if not args.sim:
			os.mkdir(maindir)
	if os.path.exists(maindir):
		os.chdir(maindir)

	# make data dir
	if not os.path.exists(DATASUBDIR):
		info('creating directory %s/%s' % (maindir, DATASUBDIR))
		if not args.sim:
			os.mkdir(DATASUBDIR)
#	elif not args.replace:
#		warning('directory %s/%s exists; skipping; use --replace to override' % (maindir, DATASUBDIR))
#		return(2)
##		continue
	if os.path.exists(DATASUBDIR):
		os.chdir(DATASUBDIR)

	if args.eigenstrat:
		opref = 'eig'
		outfile = opref + '.ind'
	else:
		opref = 'plink'
		outfile = opref + '.bed'
	if os.path.exists(outfile) and not (args.replace or args.sim):
		error('%s exists; use replace' % outfile)
		return(2)

	# make data files
	if args.s1name == 'vcf':
		if args.sepvcfs:
			for vcf, repvcf in vcf_input:
				vcfname = os.path.basename(vcf).replace('.bgz', '').replace('.vcf', '')
				jobname = ':'.join((sname, vcfname))
				if samp == '*':
					bcfview = 'bcftools view'
				else:
					bcfview = 'bcftools view -s %s ' % (samp)
				if args.eigenstrat:
					cmd = '%s %s | vcf-proc.py -H --vars | varsites2eigenstrat.py --append -p %s' % (bcfview, vcf, opref)
				else:
					cmd = '%s %s | plink --vcf /dev/stdin --double-id --allow-extra-chr --set-missing-snp-ids \"@:#\" --make-bed --out %s' % (bcfview, vcfname)
				info('processing %s' % (vcf))
				aosutils.subcall(cmd, args.sim, wait = True)
		else:
			vcfs =  ' '.join([x[0] for x in vcf_input])
			if samp == '*':
				bcfview = 'bcftools concat %s' % vcfs
			else:
				bcfview = 'bcftools concat %s | bcftools view -s %s - ' % (vcfs, samp)
			if args.eigenstrat:
				jobname = ':'.join((sname, 'vcf2eig'))
				cmd = 'bsub.py "%s | vcf-proc.py -H --vars | varsites2eigenstrat.py -p %s" -o %s.out -M 2 -j %s' % (bcfview, opref, opref, jobname)
			else:
				jobname = ':'.join((sname, 'vcf2plink'))
#				cmd = 'bsub.py "%s | set-missing-snp-ids.awk | plink --vcf /dev/stdin --double-id --allow-extra-chr --out %s" -o %s.out -M 2 -j %s' % (bcfview, opref, opref, jobname)
				cmd = 'bsub.py "%s | plink --vcf /dev/stdin --double-id --set-missing-snp-ids \"@:#\" --allow-extra-chr --make-bed --out %s" -o %s.out -M 2 -j %s' % (bcfview, opref, opref, jobname)
			if args.replace:
				cmd += ' --replace'
			info('submitting \'%s\'' % (jobname))
			aosutils.subcall(cmd, args.sim, wait = True)

	elif args.s1name == 'ms':
		jobname = ':'.join(('ms2eigenstrat', sname))
		cmd = 'ms2eigenstrat.py %s --chrlen=%d' % (ms_file, args.chrlen)
		info('running \'%s\'' % (jobname))
		aosutils.subcall(cmd, args.sim, wait = True)

def ldprune(args): # run pca
	os.chdir(args.DIR)
	debug('In %s:' % args.DIR)

	ipref = 'plink'
	opref = 'LDP_' + aosutils.fnum(args.r2)

	if os.path.exists(DATASUBDIR):
		os.chdir(DATASUBDIR)
	else:
		error('no %s dir' % DATASUBDIR)
		return(2)

	cmd = 'plink --indep-pairwise 1000kb 100 %.2f --bfile %s --out %s --allow-extra-chr' % (args.r2, ipref, opref)
	info('pruning data in %s' % (args.DIR))
	aosutils.subcall(cmd, args.sim, wait = True)
	cmd = 'plink --extract %s.prune.in --bfile %s --out %s --allow-extra-chr --make-bed' % (opref, ipref, opref)
	info('extracting non-LD SNPs in %s' % (args.DIR))
	aosutils.subcall(cmd, args.sim, wait = True)
	
def run(args): # run pca
	os.chdir(args.DIR)
	debug('In %s:' % args.DIR)

	if args.eigenstrat:
		opref = 'eig'

		if args.parfile:
			args.runpref = os.path.splitext(args.parfile)[0]

		if not args.parfile:
			args.parfile = args.runpref + '.par'
			info('writing params to %s' % args.parfile)
			if not args.sim:
				parf = open(args.parfile, 'w')
				parf.write('genotypename:\tdata/%s.geno\n' % opref)
				parf.write('snpname:\tdata/%s.snp\n' % opref)
				parf.write('indivname:\tdata/%s.ind\n' % opref)
				parf.write('evecoutname:\t%s.evec\n' % args.runpref)
				parf.write('evaloutname:\t%s.eval\n' % args.runpref)
				parf.write('numoutevec:\t%d\n' % args.numoutevec)
				parf.write('nsnpldregress:\t%d\n' % args.nsnpldregress)
				parf.write('numthreads:\t%d\n' % args.threads)
				parf.close()

		outf = args.runpref + '.out'
		jobname = ':'.join(('smartpca', args.runpref, os.path.basename(os.path.abspath(os.path.normpath(args.DIR)))))
		if not args.memory:
			args.memory = 10
		cmd = 'bsub.py "smartpca -p %s" -o %s -M %d -t %d -q %s -j %s' % (args.parfile, outf, args.memory, args.threads, args.queue, jobname)
	else:
		opref = 'plink'
		if args.ldp:
			opref = 'LDP_' + args.r2

		jobname = ':'.join(('plink-pca', os.path.basename(os.path.abspath(os.path.normpath(args.DIR)))))
		outf =  opref + '.out'
#		cmd = 'bsub.py "plink --pca --bfile %s/%s" -o %s.pca.out -M %d -t %d -q %s -j %s' % (DATASUBDIR, opref, outf, args.memory, args.threads, args.queue, jobname)
		cmd = 'plink --pca --bfile %s/%s --out %s --allow-extra-chr' % (DATASUBDIR, opref, opref)

	if args.replace:
		cmd += ' --replace'
	if args.bsim:
		cmd += ' --sim'
	if os.path.exists(outf) and not args.replace:
		warning('%s/%s exists; use --replace' % (args.DIR, outf))
	else:
		info('submitting \'%s\'' % (jobname))
		aosutils.subcall(cmd, args.sim, wait = True)

pp = argparse.ArgumentParser(add_help=False)
pp.add_argument('--replace', action='store_true', default = False, help = 'replace existing files')
pp.add_argument('--sim', action='store_true', default = False, help = 'dry run')
pp.add_argument('-v', '--verbose', action='store_true', default = False)#, help = 'dry run')
pp.add_argument('--debug', action='store_true', default = False, help=argparse.SUPPRESS)
pp.add_argument('--bsim', action='store_true', default = False, help=argparse.SUPPRESS)
pp.add_argument('--eigenstrat', action='store_true', default = False, help='use eigenstrat') 

p = argparse.ArgumentParser()
s = p.add_subparsers()#help='sub-command help')

p1 = s.add_parser('prep', help='prepare files for pca analysis')#, add_help=False)
s1 = p1.add_subparsers(dest='s1name')#help='sub-command help')
p11 = s1.add_parser('ms', parents=[pp])#, help='prep help')
p11.add_argument('MS_FILE', help = 'ms simulation file')
p11.add_argument('--chrlen', default = 50e6, help='sumulated chromosome length') 
p12 = s1.add_parser('vcf', parents=[pp])#, help='prep help')
p12.add_argument('VCF_FILE', nargs='*') 
p12.add_argument('-d', '--dirpref', help='directory name prefix') 
p12.add_argument('--sepvcfs', action='store_true', default = False, help='process vcfs separately (need to merge subsequently)') 
p12.add_argument('-s', '--samples', help='comma-separated list of sample names in VCF_FILE') 
p12.add_argument('-S', '--sample_list', help='file containing list of sample names (one per line)') 
p12.add_argument('-f', '--vcf_list', help='file containing a list of input vcfs (one per line)') 
p1.set_defaults(func=prep)

p2 = s.add_parser('run', parents=[pp], help='run smartpca')
p2.add_argument('DIR')
p2.add_argument('-p', '--parfile', default='', help='input parameter file (eigenstrat)') 
p2.add_argument('--runpref', default='smartpca', help='output file prefix') 
p2.add_argument('--numoutevec', type=int, default=4, help = 'number of output eigenvectors (smartpca)')
p2.add_argument('--nsnpldregress', type=int, default=0, help = 'number of LD regression SNPs (smartpca)')
p2.add_argument('-t', '--threads', type=int, default=8, help = 'number of threads to use')
p2.add_argument('-q', '--queue', default='normal', help = 'queue to use')
p2.add_argument('-M', '--memory', type=int, default=0, help = 'GB of RAM to use')
p2.add_argument('--ldp', action='store_true', default = False, help='use LD-pruned data') 
p2.add_argument('--r2', default = '0.2', help = 'R^2 threshold for LD-pruned data')
p2.set_defaults(func=run)

p3 = s.add_parser('ldprune', parents=[pp], help='ldprune (plink)')
p3.add_argument('DIR')
p3.add_argument('--r2', type=float, default=0.2, help = 'R^2 threshold for LD pruning')
p3.set_defaults(func=ldprune)

args = p.parse_args()

loglevel = logging.WARNING
if args.verbose:
	loglevel = logging.INFO
if args.debug:
	loglevel = logging.DEBUG
logging.basicConfig(format = '%(module)s:%(lineno)d:%(levelname)s: %(message)s', level = loglevel)

args.func(args)

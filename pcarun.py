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

EPREF = 'eig'

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
#		sname = '-'.join([os.path.splitext(os.path.basename(x))[0] for x in samp.split(',')])
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

#		# make snp-geno dir
#		subdir = 'snp'
#		if not os.path.exists(subdir):
#			info('creating directory %s/%s' % (maindir, subdir))
#			if not args.sim:
#				os.mkdir(subdir)
#		elif not args.replace:
#			warning('directory %s/%s exists; skipping; use --replace to override' % (maindir, subdir))
##			return(2)
#			continue
#		if os.path.exists(subdir):
#			os.chdir(subdir)

	# make eigenstrat files
	outfile = EPREF + '.ind'
	if os.path.exists(outfile) and not (args.replace or args.sim):
		error('%s exists; use replace' % outfile)
		return(2)

	if args.s1name == 'vcf':
		if args.vcf_list and not (args.sim or os.path.exists('tmp')):
			os.mkdir('tmp')

		for vcf, repvcf in vcf_input:
			vcfname = os.path.basename(vcf).replace('.bgz', '')#.replace('.vcf', '')
			jobname = ':'.join((sname, vcfname))
			if samp == '*':
				bcfview = 'bcftools view'
			else:
				bcfview = 'bcftools view -s %s' % (samp)
			if repvcf:
				tmprep = 'tmp/%s.rep.gz' % vcfname
				if not args.usetmp:
					cmd = '%s %s | gzip > %s' % (bcfview, repvcf, tmprep)
					info('extracting replacement calls from %s' % (repvcf))
					aosutils.subcall(cmd, args.sim, wait = True)
				reparg = '--replacecalls=%s' % tmprep
			else:
				reparg = ''
#			cmd = 'bsub.py "%s %s | vcf-proc.py -H --vars %s | varsites2eigenstrat.py --append -p %s" -o %s.out -M 2 -j %s' % (bcfview, vcf, reparg, EPREF, EPREF, jobname)
			cmd = '%s %s | vcf-proc.py -H --vars %s | varsites2eigenstrat.py --append -p %s' % (bcfview, vcf, reparg, EPREF)
			info('processing %s' % (vcf))
#			info('running \'%s\'' % (jobname))
			aosutils.subcall(cmd, args.sim, wait = True)
	elif args.s1name == 'ms':
		jobname = ':'.join(('ms2eigenstrat', sname))
		cmd = 'ms2eigenstrat.py %s --chrlen=%d' % (ms_file, args.chrlen)
		info('running \'%s\'' % (jobname))
		aosutils.subcall(cmd, args.sim, wait = True)

def run(args): # run pca
	os.chdir(args.DIR)
	debug('In %s:' % args.DIR)

	if not args.sim:
		parname = args.runpref + '.par'
		parfile = open(parname, 'w')
		parfile.write('genotypename:\t%s.geno\n' % EPREF)
		parfile.write('snpname:\t%s.snp\n' % EPREF)
		parfile.write('indivname:\t%s.ind\n' % EPREF)
		parfile.write('evecoutname:\t%s.evec\n' % args.runpref)
		parfile.write('evaloutname:\t%s.eval\n' % args.runpref)
		parfile.write('numoutevec:\t%d\n' % args.numoutevec)
		parfile.write('nsnpldregress:\t%d\n' % args.nsnpldregress)
		parfile.write('numthreads:\t%d\n' % args.threads)
		parfile.close()

	outf = args.runpref + '.out'
	jobname = ':'.join(('smartpca', args.runpref, os.path.basename(os.path.abspath(os.path.normpath(args.DIR)))))
	if not args.memory:
		args.memory = 10
	cmd = 'bsub.py "smartpca -p %s" -o %s -M %d -t %d -q %s -j %s' % (parname, outf, args.memory, args.threads, args.queue, jobname)
	if args.replace:
		cmd += ' --replace'
	if args.bsim:
		cmd += ' --sim'
	if os.path.exists(outf) and not args.replace:
		warning('%s/%s exists; use --replace' % (args.DIR, outf))
	else:
		info('submitting \'%s\'' % (jobname))
		aosutils.subcall(cmd, args.sim, wait = True)

def plot(args): # make plots
	if not args.sim:
		os.chdir(args.DIR)
	sname = os.path.basename(os.path.normpath(args.DIR))
	jobname = ':'.join(('smcplot', sname))

	cmd = 'smcplot.py -m run.final.txt -t %f -u %e --hc -o %s' % (args.tgen, args.mu, sname)
	info('running \'%s\'' % (jobname))
	aosutils.subcall(cmd, args.sim, wait = True)


pp = argparse.ArgumentParser(add_help=False)
pp.add_argument('--replace', action='store_true', default = False, help = 'replace existing files')
pp.add_argument('--sim', action='store_true', default = False, help = 'dry run')
pp.add_argument('-v', '--verbose', action='store_true', default = False)#, help = 'dry run')
pp.add_argument('--debug', action='store_true', default = False, help=argparse.SUPPRESS)
pp.add_argument('--bsim', action='store_true', default = False, help=argparse.SUPPRESS)

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
p12.add_argument('-s', '--samples', help='comma-separated list of sample names in VCF_FILE') 
p12.add_argument('-S', '--sample_list', help='file containing list of sample names (one per line)') 
p12.add_argument('-f', '--vcf_list', help='file containing a list of input vcfs (one per line). For each one, a vcf of replacement calls may specified in a second column.') 
p12.add_argument('--usetmp', action='store_true', default=False, help='use existing replacement calls in tmp dir') 
p1.set_defaults(func=prep)

p2 = s.add_parser('run', parents=[pp], help='run smartpca')
p2.add_argument('DIR')
p2.add_argument('-p', '--runpref', default='eig', help='output file prefix') 
p2.add_argument('--numoutevec', type=int, default=4, help = 'number of output eigenvectors (smartpca)')
p2.add_argument('--nsnpldregress', type=int, default=0, help = 'number of LD regression SNPs (smartpca)')
p2.add_argument('-t', '--threads', type=int, default=8, help = 'number of threads to use')
p2.add_argument('-q', '--queue', default='normal', help = 'queue to use')
p2.add_argument('-M', '--memory', type=int, default=0, help = 'GB of RAM to use')
p2.set_defaults(func=run)

#p3 = s.add_parser('plot', parents=[pp], help='plot pca results')
#p3.add_argument('DIR')
##p3.add_argument('-t', '--tgen', type=int, default=30.0, help = 'generation time (y)')
#p3.set_defaults(func=plot)

args = p.parse_args()

loglevel = logging.WARNING
if args.verbose:
	loglevel = logging.INFO
if args.debug:
	loglevel = logging.DEBUG
logging.basicConfig(format = '%(module)s:%(lineno)d:%(levelname)s: %(message)s', level = loglevel)

args.func(args)

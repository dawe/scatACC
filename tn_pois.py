import sys
import pysam
import pyBigWig
import numpy as np
import scipy.stats as sst
import argparse

_MAXLOG = 64 #reasonable value?

def get_options():
	parser = argparse.ArgumentParser(prog='tn_pois.py')
	parser.add_argument('-p', '--prefix', help='Prefix for output files', default='out')
	parser.add_argument('-s', '--stepsize', help='Bin size for output', default=5000, type=int)
	parser.add_argument('-5', '--tn5', help='BAM file with tn5 signal', required=True)
	parser.add_argument('-H', '--tnH', help='BAM file with tnH signal', required=True)
	parser.add_argument('-A', '--accurate', help='Accurate coverage estimate (slow)', action='store_true')
	parser.add_argument('-M', '--mito_chrom_name', help='Name of mitochondrial chromosome to be excluded', default='chrM')
	
	options = parser.parse_args()
	
	return options

def genome_size(fh, mt_name='chrM'):
	mt_len = fh.get_reference_length(mt_name)
	g_len = sum([fh.get_reference_length(r) for r in fh.references])
	return g_len - mt_len

def mapped_reads(fh, mt_name='chrM'):
	# calculate genome size
	references = fh.references
	nref = len(references)
	stats = fh.get_index_statistics()
	return sum([stats[x].mapped for x in range(nref) if references[x] != mt_name])


def diff_estimate():
	options = get_options()
	if options.accurate:
		# to be moved downward sometime
		raise NotImplementedError("Accurate estimate not yet implemented!")
		sys.exit(1)
	
	ftn5 = pysam.AlignmentFile(options.tn5, 'rb')
	ftnH = pysam.AlignmentFile(options.tnH, 'rb')
	
	n_reads5 = mapped_reads(ftn5, options.mito_chrom_name)
	n_readsH = mapped_reads(ftnH, options.mito_chrom_name)
	
	# genome sizes should be the same, but you'll never know
	# at a certain point one should also check that two files
	# have been aligned on the same reference...
	gs5 = genome_size(ftn5, options.mito_chrom_name)
	gsH = genome_size(ftnH, options.mito_chrom_name)
	
	# poisson lambda will be the average coverage per bin
	mu5 = n_reads5 / (gs5 / options.stepsize)
	muH = n_readsH / (gsH / options.stepsize)	
	
	# now we can use a Skellam distribution
	# this will model the differences we count on each bin
	
	skd = sst.skellam(mu5, muH)
	
	# by default we will skip duplicated reads
	for chromosome in ftn5.references:
		if chromosome == options.mito_chrom_name:
			continue
		chr_len = ftn5.get_reference_length(chromosome)
		for start in range(0, chr_len, options.stepsize):
			stop = start + options.stepsize
			if stop > chr_len:
				stop = chr_len
			c5 = ftn5.count(contig=chromosome, start=start, stop=stop, read_callback='all')
			cH = ftnH.count(contig=chromosome, start=start, stop=stop, read_callback='all')	
			d = c5 - cH
			p = skd.cdf(d)
			l = np.log(p)
			if p > 0.5:
				# more evidences on the right tail, the tn5 one
				l = -np.log(1 - p)
			if l > _MAXLOG:
				l = _MAXLOG
			if l < -_MAXLOG:
				l = -_MAXLOG
			# write to stdout
			# one day I will add direct bigwig support
			sys.stdout.write(f'{chromosome}\t{start}\t{stop}\t{l:.5e}\n')
			
	
		

if __name__ == '__main__':
	diff_estimate()

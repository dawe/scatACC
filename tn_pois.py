import sys
import pysam
import pyBigWig
import scipy.stats as sst
import argparse

def get_options():
	parser = argparse.ArgumentParser(prog='tn_pois.py')
	parser.add_argument('-p', '--prefix', help='Prefix for output files', default='out')
	parser.add_argument('-s', '--stepsize', help='Bin size for output', default=5000, type=int)
	parser.add_argument('-5', '--tn5', help='BAM file with tn5 signal', required=True)
	parser.add_argument('-H', '--tnH', help='BAM file with tnH signal', required=True)
	parser.add_argument('-A', '--accurate', help='Accurate coverage estimate (slow)', action='store_true')	
	
	options = parser.parse_args()
	
	return options



def diff_estimate():
	options = get_options()
	if options.accurate:
		raise NotImplementedError("Accurate estimate not yet implemented!")
		sys.exit(1)
	

if __name__ == '__main__':
	diff_estimate()

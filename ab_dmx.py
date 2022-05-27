import HTSeq
import sys
import argparse
import editdistance as ed
import json

_SPACER_MAXDIST = 2
_BC_MAXDIST = 1


bc_A = ['CGTACTAG','TCCTGAGC', 'TCATGAGC', 'CCTGAGAT', 
'TAAGGCGA', 'GCTACGCT', 'AGGCTCCG', 'CTGCGCAT']
bc_B = ['AACC', 'GCTC']

_bc_A = [b'CGTACTAG',b'TCCTGAGC', b'TCATGAGC', b'CCTGAGAT', 
b'TAAGGCGA', b'GCTACGCT', b'AGGCTCCG', b'CTGCGCAT']
_bc_B = [b'AACC', b'GCTC']


_sp_A = b'AGATGTGTATAAGAGACAG'

dbc_A = {
'CGTACTAG':'tn5',
'TCCTGAGC':'tn5',
'TCATGAGC':'tn5',
'CCTGAGAT':'tn5',
'TAAGGCGA':'tnH',
'GCTACGCT':'tnH',
'AGGCTCCG':'tnH',
'CTGCGCAT':'tnH'
}

_dbc_A = {
b'CGTACTAG':'CGTACTAG',
b'TCCTGAGC':'TCCTGAGC',
b'TCATGAGC':'TCATGAGC',
b'CCTGAGAT':'CCTGAGAT',
b'TAAGGCGA':'TAAGGCGA',
b'GCTACGCT':'GCTACGCT',
b'AGGCTCCG':'AGGCTCCG',
b'CTGCGCAT':'CTGCGCAT'
}


dbc_B = {
'AACC':'tn5',
'GCTC':'tnH'
}

_dbc_B = {
b'AACC':'AACC',
b'GCTC':'GCTC'
}


def get_options():
	parser = argparse.ArgumentParser(prog='ab_tmx.py')
	parser.add_argument('-p', '--prefix', help='Prefix for output files', default='out')
	parser.add_argument('-1', '--read1', help='Fastq file with R1', required=True)
	parser.add_argument('-2', '--read2', help='Fastq file with R2 (is R3 in scGETseq)', required=True)
#	parser.add_argument('-d', '--definitions', help='JSON file with barcode definitions', required=True)	
	parser.add_argument('-b', '--barcodes', help='Fastq file with cell barcodes (is R2 in scGETseq only)')	
	parser.add_argument('-U', '--write_unmatched', help='Dump unmatched reads', store_action=True)	
	
	options = parser.parse_args()
	
	return options

def demux():
	options = get_options()	
	out_r1 = 'R1'
	out_r2 = 'R2'
	suffix = [out_r1, out_r2]
	if options.barcodes:
		suffix = [out_r1, out_r2, 'RB']
		
	filenames = [f'{options.prefix}_{a}_{b}_{r}_{dbc_A[a]}_{dbc_B[b]}.fastq' for a in bc_A for b in bc_B for r in suffix]

	wfh = dict.fromkeys(filenames)
	for f in wfh:
		wfh[f] = open(f, 'w')
	
	r1 = HTSeq.FastqReader(options.read1)
	r2 = HTSeq.FastqReader(options.read2)	
	if options.barcodes:
		rb = HTSeq.FastqReader(options.barcodes)
		read_iterator = zip(r1, r2, rb)
	else:
		read_iterator = zip(r1, r2)
	
	n_tot = 0
	n_pass = 0
	for reads in read_iterator:
		# check MEDSA spacer
		n_tot += 1
		if ed.eval(reads[0].seq[8:27], _sp_A) > _SPACER_MAXDIST:
			continue
		bc1 = [x for x in _bc_A if ed.eval(reads[0].seq[:8], x) <= 1]
		bc2 = [x for x in _bc_B if ed.eval(reads[1].seq[:4], x) <= 1]
		
		if len(bc1) == 1 and len(bc2) == 1:
			n_pass += 1
			bc1 = _dbc_A[bc1[0]]
			bc2 = _dbc_B[bc2[0]]
			fname1 = f'{options.prefix}_{bc1}_{bc2}_R1_{dbc_A[bc1]}_{dbc_B[bc2]}.fastq'
			fname2 = f'{options.prefix}_{bc1}_{bc2}_R2_{dbc_A[bc1]}_{dbc_B[bc2]}.fastq'
			# found
			reads[0][27:].write_to_fastq_file(wfh[fname1])
			reads[1][4:].write_to_fastq_file(wfh[fname2])
			if options.barcodes:
				fnameb = f'{options.prefix}_{bc1}_{bc2}_RB_{dbc_A[bc1]}_{dbc_B[bc2]}.fastq'
				reads[2].write_to_fastq_file(wfh[fnameb])

	for f in wfh:
		wfh[f].close()
	
	eff = n_pass / n_tot * 100
	sys.stderr.write(f'Found {n_pass} out of {n_tot} sequences {eff:.3f}%\n')
	
	
	
if __name__ == '__main__':
	demux()


import HTSeq
import sys
import argparse
import editdistance as ed
import json

_SPACER_MAXDIST = 2
_BC_MAXDIST = 1


bc_A = ['CGTACTAG','TCCTGAGC', 'TCATGAGC', 'CCTGAGAT', 
'TAAGGCGA', 'GCTACGCT', 'AGGCTCCG', 'CTGCGCAT']
bc_B = ['AACCAACC', 'GCTCGCTC']

_bc_A = [x.encode() for x in bc_A]
_bc_B = [x.encode() for x in bc_B]

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

dbc_B = {
'AACCAACC':'tn5',
'GCTCGCTC':'tnH'
}

_dbc_A = dict(zip(_bc_A, bc_A))
_dbc_B = dict(zip(_bc_B, bc_B))


def get_options():
	parser = argparse.ArgumentParser(prog='ab_tmx.py')
	parser.add_argument('-p', '--prefix', help='Prefix for output files', default='out')
	parser.add_argument('-1', '--read1', help='Fastq file with R1', required=True)
	parser.add_argument('-2', '--read2', help='Fastq file with R2 (is R3 in scGETseq)', required=True)
#	parser.add_argument('-d', '--definitions', help='JSON file with barcode definitions', required=True)	
	parser.add_argument('-b', '--barcodes', help='Fastq file with cell barcodes (is R2 in scGETseq only)')	
	parser.add_argument('-U', '--write_unmatched', help='Dump unmatched reads', action='store_true')	
	
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
	if options.write_unmatched:
		filenames = filenames + [f'{options.prefix}_un_{r}.fastq' for r in suffix]

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
		# check MEDSA  and MEDSB spacer
		n_tot += 1
		if ed.eval(reads[0].seq[8:27], _sp_A) > _SPACER_MAXDIST and ed.eval(reads[1].seq[8:27], _sp_A) > _SPACER_MAXDIST:
			if options.write_unmatched:
				fname1 = f'{options.prefix}_un_R1.fastq'
				fname2 = f'{options.prefix}_un_R2.fastq'
				reads[0].write_to_fastq_file(wfh[fname1])
				reads[1].write_to_fastq_file(wfh[fname2])
				if options.barcodes:
					fnameb = f'{options.prefix}_un_RB.fastq'
					reads[2].write_to_fastq_file(wfh[fnameb])
			continue
		bc_dist1 = [ed.eval(reads[0].seq[:8], x)for x in _bc_A]
		amin1 = [x for x in range(len(bc_dist1)) if bc_dist1[x] == min(bc_dist1)][0]
		bc_dist2 = [ed.eval(reads[1].seq[:8], x)for x in _bc_B]
		amin2 = [x for x in range(len(bc_dist2)) if bc_dist2[x] == min(bc_dist2)][0]
	
		if bc_dist1[amin1] <= _BC_MAXDIST and bc_dist2[amin2] <= _BC_MAXDIST:
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


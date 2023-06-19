import HTSeq
import sys
import argparse
#import editdistance as ed
import bgzip


_SPACER_MAXDIST = 2
_BC_MAXDIST = 1

## GETseq barcode definitions
bc_A = ['CGTACTAG','TCCTGAGC', 'TCATGAGC', 'CCTGAGAT', 
'TAAGGCGA', 'GCTACGCT', 'AGGCTCCG', 'CTGCGCAT']

# binary version to make things quicker
_bc_A = [x.encode() for x in bc_A]

#MEDS spacer
_sp_A = b'AGATGTGTATAAGAGACAG'

# tn association
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

# binary decoder to be quicker
_dbc_A = dict(zip(_bc_A, bc_A))

def hamming(a, b):
    la = len(a)
    lb = len(b)
    if la != lb:
        return la
    return sum([a[x] != b[x] for x in range(len(a))])

def empty_spool(d):
    for k in d:
        d[k] = b''

def get_options():
    parser = argparse.ArgumentParser(prog='ab_tmx.py')
    parser.add_argument('-p', '--prefix', help='Prefix for output files', default='out')
    parser.add_argument('-1', '--read1', help='Fastq file with R1', required=True)
    parser.add_argument('-2', '--read2', help='Fastq file with R2 (is R3 in scGETseq)', required=True)
    parser.add_argument('-b', '--barcodes', help='Fastq file with cell barcodes (is R2 in scGETseq only)')    
    parser.add_argument('-U', '--write_unmatched', help='Dump unmatched reads', action='store_true')    
    parser.add_argument('-T', '--tagdust', help='tagdust compatible naming (for retrocompatibility with older getseq pipeline)', action='store_true')
    parser.add_argument('-n', '--n_seq', help='Max number of sequences to process (for debugging)', default=0, type=int)

    
    options = parser.parse_args()
    
    return options

def demux():
    nl = b'\n'
    dnl = b'\n+\n'
    _chunk_size = 512 # number 

    options = get_options()    
    out_r1 = 'R1'
    out_r2 = 'R2'
    suffix = [out_r1, out_r2]
    if options.barcodes:
        suffix = [out_r1, out_r2, 'RB']

    filenames = [f'{options.prefix}_{a}_{r}_{dbc_A[a]}.fastq.gz' for a in bc_A for r in suffix]
    un_filenames = []
    if options.write_unmatched:
        un_filenames =  [f'{options.prefix}_un_{r}.fastq.gz' for r in suffix]
        
    if options.tagdust:
        suffix = ['READ1', 'READ3', 'READ2']
        filenames = [f'{options.prefix}_BC_{a}_{r}.fq.gz' for a in bc_A for r in suffix]        
        un_filenames = [f'{options.prefix}_un_{r}.fq.gz' for r in suffix]
        options.write_unmatched = True

    wfh = dict.fromkeys(filenames + un_filenames)
    wfh_bgz = dict.fromkeys(filenames + un_filenames)
    spool = dict.fromkeys(filenames + un_filenames)
    for f in wfh:
        wfh[f] = open(f, 'wb')
        wfh_bgz[f] = bgzip.BGZipWriter(wfh[f])
        spool[f] = b''
    
    r1 = HTSeq.FastqReader(options.read1)
    r2 = HTSeq.FastqReader(options.read2)    
    if options.barcodes:
        rb = HTSeq.FastqReader(options.barcodes)
        read_iterator = zip(r1, r2, rb)
    else:
        read_iterator = zip(r1, r2)
    
    n_tot = 0
    n_pass = 0
    _spool_counter = 0
    for reads in read_iterator:
        # check MEDSA spacer
        
        if options.n_seq > 0 and n_tot == options.n_seq:
            break
        
        seq1 = reads[0].seq
        seq2 = reads[1].seq
        seqb = reads[2].seq
        
        qual1 = reads[0].qualstr
        qual2 = reads[1].qualstr
        qualb = reads[2].qualstr
        
        name1 = f'@{reads[0].name}'.encode()
        name2 = f'@{reads[1].name}'.encode()
        nameb = f'@{reads[2].name}'.encode()
        
        n_tot += 1
        if hamming(seq1[8:27], _sp_A) > _SPACER_MAXDIST:
            if options.write_unmatched:
                spool[un_filenames[0]] = spool[un_filenames[0]] + name1 + nl + seq1 + dnl + qual1 + nl
                spool[un_filenames[1]] = spool[un_filenames[1]] + name2 + nl + seq2 + dnl + qual2 + nl
                if options.barcodes:
                    spool[un_filenames[2]] = spool[un_filenames[2]] + nameb + nl + seqb + dnl + qualb + nl
            continue    
        bc_dist = [hamming(seq1[:8], x)for x in _bc_A]
        amin = [x for x in range(len(bc_dist)) if bc_dist[x] == min(bc_dist)][0]
        
        if bc_dist[amin] <= _BC_MAXDIST:
            n_pass += 1
            bc1 = _dbc_A[_bc_A[amin]]
            if options.tagdust:
                fname1 = f'{options.prefix}_BC_{bc1}_READ1.fq.gz'
                fname2 = f'{options.prefix}_BC_{bc1}_READ3.fq.gz'
            else:
                fname1 = f'{options.prefix}_{bc1}_R1_{dbc_A[bc1]}.fastq.gz'
                fname2 = f'{options.prefix}_{bc1}_R2_{dbc_A[bc1]}.fastqgz'
            # found
            spool[fname1] = spool[fname1] + name1 + nl + seq1[27:] + dnl + qual1[27:] + nl
            spool[fname2] = spool[fname2] + name2 + nl + seq2 + dnl + qual2 + nl
            if options.barcodes:
                if options.tagdust:
                    fnameb = f'{options.prefix}_BC_{bc1}_READ2.fq.gz'
                else:
                    fnameb = f'{options.prefix}_{bc1}_RB_{dbc_A[bc1]}.fastq.gz'
                spool[fnameb] = spool[fnameb] + nameb + nl + seqb + dnl + qualb + nl    

        elif options.write_unmatched:        
            spool[un_filenames[0]] = spool[un_filenames[0]] + name1 + nl + seq1 + dnl + qual1 + nl
            spool[un_filenames[1]] = spool[un_filenames[1]] + name2 + nl + seq2 + dnl + qual2 + nl
            if options.barcodes:
                spool[un_filenames[2]] = spool[un_filenames[2]] + nameb + nl + seqb + dnl + qualb + nl

        _spool_counter += 1
        if _spool_counter == _chunk_size:
            for k in wfh_bgz.keys():
                wfh_bgz[k].write(spool[k])
                spool[k] = b''
            _spool_counter = 0

    # end, write remaining spool and close files
    for k in wfh_bgz.keys():
        if len(spool[k]) > 0:
            wfh_bgz[k].write(spool[k])
            wfh_bgz[k].close()
        wfh[k].close()
        

    eff = n_pass / n_tot * 100
    sys.stderr.write(f'Found {n_pass} out of {n_tot} sequences {eff:.3f}%\n')
    
    
    
if __name__ == '__main__':
    demux()


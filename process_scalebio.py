import pandas as pd
import sys
import numpy as np
import bgzip
#import editdistance as ed
from tqdm import tqdm
import argparse
import HTSeq

# LLLLLLLLL[L] CAGAGC UUUUUUUU RRRRRRRRRR polyT
# L: ligation barcode (can be 9 or 10 nt long)
# U: Umi
# R: RT barcode

def hamming(a, b):
    la = len(a)
    lb = len(b)
    if la != lb:
        return la
    return sum([a[x] != b[x] for x in range(len(a))])

def correct_bc(bc, bc_list, distance=1):
    distances = np.array([hamming(bc, sh_bc) for sh_bc in bc_list])
    accepted = np.where(distances <= distance)[0]
    if len(accepted) >0:
        # return the first one
        return bc_list[accepted[0]]
    else:
        return b''
 
def get_options():
    parser = argparse.ArgumentParser(prog='process_share.py')
    parser.add_argument('-1', '--read1', help='Read 1, containing cell barcode (R1)', required=True)
    parser.add_argument('-2', '--read2', help='Read 2 (R2)', required=True)
    parser.add_argument('-F', '--filter_failed', help='Filter failed reads (without proper structure)', action='store_true')
    parser.add_argument('-p', '--prefix', help='Prefix for output files')
    parser.add_argument('-R', '--rt_correct_file', help='List of RT barcodes', default='', required=True)
    parser.add_argument('-L', '--lig_correct_file', help='List of ligation barcodes', default='', required=True)
    parser.add_argument('-t', '--threshold', help='Max distance when correcting barcodes', default=1, type=int)
    parser.add_argument('-n', '--n_seq', help='Max number of sequences to process (for debugging)', default=0, type=int)
  
    options = parser.parse_args()
  
    return options

def main():
    nl = b'\n'
    dnl = b'\n+\n'
    linker = b'CAGAGC'
    _chunk_size = 512 # number 

    options = get_options()
    
    bc_fix = True
    rt_list = []
    lig_list = []
    for line in open(options.rt_correct_file):
        t = line.split()
        rt_list.append(bytes(t[0], encoding='ascii'))
    rt_list = np.array(rt_list)

    for line in open(options.lig_correct_file):
        t = line.split()
        lig_list.append(bytes(t[0], encoding='ascii'))
    lig_list = np.array(lig_list)
    
    
    sp1 = 'CAGAGC' # [9:15] or [10:16]

    r1 = HTSeq.FastqReader(options.read1)
    r2 = HTSeq.FastqReader(options.read2)

    read_iterator = zip(r1, r2)
    
    fout1 = f'{options.prefix}_R1.fastq.gz'
    fout2 = f'{options.prefix}_R2.fastq.gz'
    
    raw_out1 = open(fout1, 'wb')
    raw_out2 = open(fout2, 'wb')
    fh_out1 = bgzip.BGZipWriter(raw_out1, batch_size=256)
    fh_out2 = bgzip.BGZipWriter(raw_out2, batch_size=256)
    # remember to write bytes, not strings

    r1_spool = b''
    r2_spool = b''
    _spool_counter = 0
    
    n_tot = 0
    n_fail = 0
    for item in read_iterator:
        if options.n_seq > 0 and n_tot == options.n_seq:
            break
    
        seq1 = item[0].seq
        seq2 = item[1].seq

        qual1 = item[0].qualstr
        qual2 = item[1].qualstr

        n_tot += 1
       
        if options.filter_failed and not (seq1[9:15] == linker or seq1[10:16] == linker):
            n_fail += 1
            continue

        name1 = bytes('@' + item[0].name, encoding='ascii')
        name2 = bytes('@' + item[1].name, encoding='ascii')

        bc1 = seq1[:9]
        u_start = 15
        if seq1[9:15] == linker:
            u_start = 15
        elif seq1[10:16] == linker:
            u_start = 16
        umi = seq1[u_start:(u_start + 8)]
        bc2 = seq1[(u_start + 8):(u_start + 18)]
        
        q_bc1 = qual1[:9]
        q_umi = qual1[u_start:(u_start + 8)]
        q_bc2 = qual1[(u_start + 8):(u_start + 18)]
        
        if bc_fix:
            bc1 = correct_bc(bc1, lig_list, options.threshold)
            bc2 = correct_bc(bc2, rt_list, options.threshold)
            
        seq1_out = bc1 + bc2 + umi
        q_seq1_out = q_bc1 + q_bc2 + q_umi
        
        # since bc correction returns empty bc if not found
        # we can use it to skip bad reads
        
        if options.filter_failed and len(seq1_out) != len(q_seq1_out):
            n_fail += 1
            continue
        
        r1_spool = r1_spool + name1 + nl + seq1_out + dnl + q_seq1_out + nl
        r2_spool = r2_spool + name2 + nl + seq2 + dnl + qual2 + nl

        _spool_counter += 1
        
        if _spool_counter == _chunk_size:
            fh_out1.write(r1_spool)
            fh_out2.write(r2_spool)
            r1_spool = b''
            r2_spool = b''
            _spool_counter = 0
    
    # end, write remaining spool and close files
    fh_out1.write(r1_spool)
    fh_out1.close()
    raw_out1.close()
    fh_out2.write(r2_spool)
    fh_out2.close()
    raw_out2.close()

    n_pass = n_tot - n_fail
    sys.stderr.write(f"Total sequences:\t{n_tot}\n")
    f = n_fail / n_tot * 100
    sys.stderr.write(f"Failed BC:\t{n_fail} ({f:.3f}%)\n")
    f = n_pass / n_tot * 100
    sys.stderr.write(f"Passing sequences:\t{n_pass} ({f:.3f}%)\n")
    
#    eff = n_pass / n_tot * 100
#    sys.stderr.write(f'Found {n_pass} out of {n_tot} sequences {eff:.3f}%\n')
#    eff = n_dark / n_tot * 100
#    sys.stderr.write(f'Found {n_dark} dark sequences {eff:.3f}%\n')
#    eff = n_fail / (n_tot - n_dark) * 100
#    sys.stderr.write(f'Could not fix barcode for {n_fail} sequences {eff:.3f}%\n')



if __name__ == '__main__':
    main()

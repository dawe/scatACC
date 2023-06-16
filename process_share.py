import pandas as pd
import sys
import numpy as np
import bgzip
#import editdistance as ed
from tqdm import tqdm
import argparse
import HTSeq


# struct Barcode
# read     CCGAGCCCACGAGAC | TCGGACGATCATGGGNNNNNNNNCAAGTATGCAGCGCGCTCAAGCACGTGGATNNNNNNNNAGTCGTACGCCGATGCGAAACATCGGCCACNNNNNNNN | ATCTCGTATGCCGTCTTCTGCTTG
# template GGCTCGGGTGCTCTG | AGCCTGCTAGTACCCNNNNNNNNGTTCATACGTCGCGCGAGTTCGTGCACCTANNNNNNNNTCAGCATGCGGCTACGCTTTGTAGCCGGTGNNNNNNNN | TAGAGCATACGGCAGAAGACGAAC
# struct UMI
# AGATGTGTATAAGAGACAG | NNNNNNNNNN


def get_options():
    parser = argparse.ArgumentParser(prog='process_share.py')
    parser.add_argument('-1', '--read1', help='Read 1 (R1)')
    parser.add_argument('-2', '--read2', help='Read 2, containing cell barcode (R2)')
    parser.add_argument('-3', '--read3', help='Read 3, containing either UMI or second ATAC pair (R3)')
    parser.add_argument('-R', '--rna', help='Process as scRNA-seq, stitching R3 and R2', action='store_true')
    parser.add_argument('-L', '--stitch_length', help='Number of bp to retain when stitching', default=10)
    parser.add_argument('-F', '--filter_failed', help='Filter failed reads (having GGGGG stretches)', action='store_true')
    parser.add_argument('-p', '--prefix', help='Prefix for output files')
    parser.add_argument('-C', '--bc_correct_file', help='Fix cell barcode to given list	')
    parser.add_argument('-t', '--threshold', help='Max distance when correcting barcodes', default=1, type=int)
  
    options = parser.parse_args()
  
    return options

def main():
    nl = b'\n'
    dark = b'GGGGGGGGGGGGGGGGGGGG'
    _chunk_size = 2048 # number 

    options = get_options()
    
    sp1 = 'TCGGACGATCATGGG' # [0:15]
    sp2 = 'CAAGTATGCAGCGCGCTCAAGCACGTGGAT' # [23:53]
#    sp3 = 'AGTCGTACGCCGATGCGAAACATCGGCCAC' # [61:91]
#            AGTCGTACGCCGATGCGAAACATCGGCCAC
    sp3 = 'AGTCGTACGCCGATGCAAAACATAAACCAC' # [61:91]

    r1_out = open

    r1 = HTSeq.FastqReader(options.read1)
    r2 = HTSeq.FastqReader(options.read2)
    r3 = HTSeq.FastqReader(options.read3)

    read_iterator = zip(r1, r2, r3)
    
    fout1 = f'{options.prefix}_R1.fastq.gz'
    fout2 = f'{options.prefix}_R2.fastq.gz'
    fout3 = f'{options.prefix}_R3.fastq.gz'
    
    raw_out1 = open(fout1, 'wb')
    raw_out2 = open(fout2, 'wb')
    fh_out1 = bgzip.BGZipWriter(raw_out1)
    fh_out2 = bgzip.BGZipWriter(raw_out2)
    if not options.rna:
        raw_out3 = open(fout3, 'wb')
        fh_out3 = bgzip.BGZipWriter(raw_out3)
    # remember to write bytes, not strings

    for item in read_iterator:
        seq1 = item[0].seq
        seq2 = item[1].seq
        seq3 = item[2].seq

        qual1 = item[0].qualstr
        qual2 = item[1].qualstr
        qual3 = item[2].qualstr

        name1 = item[0].name
        name2 = item[1].name
        name3 = item[2].name
        
        is_dark = False
        if seq3[:20] == dark:
            is_dark = True
        
        if options.filter_failed and is_dark:
            continue
            
        
        
        
        
    
    if options.read_umi:
        umi_fh = iter(HTSeq.FastqReader(options.read_umi))
        _iter = zip(bc_fh, umi_fh)
    else:
        _iter = bc_fh

    for _r in _iter:
        if options.read_umi:
            bc_r, umi_r = _r
        else:
            bc_r = _r
    
        bc_s = bc_r.seq.decode('ascii')
        bc_q = bc_r.qualstr.decode('ascii')
        e1 = ed.eval(bc_s[:15], sp1)
        e2 = ed.eval(bc_s[23:53], sp2)
        e3 = ed.eval(bc_s[61:91], sp3)
        if options.debug:
            sys.stderr.write(f'{e1}\t{e2}\t{e3}\n')
        if not options.crop:
            if e1 > thr[0] and e2 > thr[1] and e3 > thr[2]:
                continue
            comb_bc = bc_s[15:23] + bc_s[53:61] + bc_s[91:99]
            comb_ql = bc_q[15:23] + bc_q[53:61] + bc_q[91:99]
            if 'N' in comb_bc:
                continue
        else:
            comb_bc = bc_s[15:23] + bc_s[53:61] + bc_s[91:99]
            comb_ql = bc_q[15:23] + bc_q[53:61] + bc_q[91:99]
        
        if options.read_umi:
            umi_s = umi_r.seq.decode('ascii')
            umi_q = umi_r.qualstr.decode('ascii')
            comb_bc = comb_bc + umi_s[:10]
            comb_ql = comb_ql + umi_q[:10]
        
        sys.stdout.write(f'@{bc_r.name}\n{comb_bc}\n+\n{comb_ql}\n')
            


if __name__ == '__main__':
    main()

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
    parser.add_argument('-1', '--read1', help='Read 1 (R1)')
    parser.add_argument('-2', '--read2', help='Read 2, containing cell barcode (R2)')
    parser.add_argument('-3', '--read3', help='Read 3, containing either UMI or second ATAC pair (R3)')
    parser.add_argument('-R', '--rna', help='Process as scRNA-seq, stitching R3 and R2', action='store_true')
    parser.add_argument('-L', '--stitch_length', help='Number of bp to retain when stitching', default=10)
    parser.add_argument('-F', '--filter_failed', help='Filter failed reads (having GGGGG stretches)', action='store_true')
    parser.add_argument('-p', '--prefix', help='Prefix for output files')
    parser.add_argument('-C', '--bc_correct_file', help='Fix cell barcode to given list	', default='')
    parser.add_argument('-t', '--threshold', help='Max distance when correcting barcodes', default=1, type=int)
  
    options = parser.parse_args()
  
    return options

def main():
    nl = b'\n'
    dnl = b'\n+\n'
    dark = b'GGGGGGGGGGGGGGGGGGGG'
    _chunk_size = 2048 # number 

    options = get_options()
    
    bc_fix = False
    if options.bc_correct_file:
        bc_fix = True
        bc_list = []
        for line in open(options.bc_correct_file):
            t = line.split()
            bc_list.append(bytes(t[1], encoding='ascii'))
        bc_list = np.array(bc_list)
    
    
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

    r1_spool = b''
    r2_spool = b''
    r3_spool = b''
    _spool_counter = 0
    
    for item in read_iterator:
        seq1 = item[0].seq
        seq2 = item[1].seq
        seq3 = item[2].seq

        qual1 = item[0].qualstr
        qual2 = item[1].qualstr
        qual3 = item[2].qualstr

        
        is_dark = False
        if seq3[:20] == dark:
            is_dark = True
        
        if options.filter_failed and is_dark:
            continue

        name1 = bytes('@' + item[0].name, encoding='ascii')
        name2 = bytes('@' + item[1].name, encoding='ascii')
        name3 = bytes('@' + item[2].name, encoding='ascii')
            
        bc1 = seq2[15:23]
        bc2 = seq2[53:61]
        bc3 = seq2[91:99]
        
        q_bc1 = qual2[15:23]
        q_bc2 = qual2[53:61]
        q_bc3 = qual2[91:99]
        
        if bc_fix:
            bc1 = correct_bc(bc1, bc_list, options.threshold)
            bc2 = correct_bc(bc2, bc_list, options.threshold)
            bc3 = correct_bc(bc3, bc_list, options.threshold)
            
        seq2_out = bc1 + bc2 + bc3
        q_seq2_out = q_bc1 + q_bc2 + q_bc3
        
        # since bccorrection returns empty bc if not found
        # we can use it to skip bad reads
        
        if len(seq2_out) != q_seq2_out and options.filter_failed:
            continue
        
        if options.rna:
            seq2_out = seq2_out + seq3[:options.stitch_length]
            q_seq2_out = q_seq2_out + qual3[:options.stitch_length]
            
        r1_spool = r1_spool + name1 + nl + seq1 + dnl + qual1 + nl
        r2_spool = r2_spool + name2 + nl + seq2_out + dnl + q_seq2_out + nl
        
        if not options.rna:
            r3_spool = r3_spool + name3 + nl + seq3 + dnl + qual3 + nl
        
        _spool_counter += 1
        
        if _spool_counter == _chunk_size:
            fh_out1.write(r1_spool)
            fh_out2.write(r2_spool)
            if not options.rna:
                fh_out3.write(r3_spool)
            r1_spool = b''
            r2_spool = b''
            r3_spool = b''
            _spool_counter = 0
    
    # end, write remaining spool and close files
    if len(r1_spool) > 0:
        fh_out1.write(r1_spool)
        fh_out2.write(r2_spool)
        if not options.rna:
            fh_out3.write(r3_spool)
            
    fh_out1.close()
    fh_out2.close()
    if not options.rna:
        fh_out3.close()


            
        
            


if __name__ == '__main__':
    main()

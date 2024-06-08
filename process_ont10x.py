import pandas as pd
import sys
import numpy as np
import bgzip
#import editdistance as ed
from tqdm import tqdm
import argparse
import HTSeq

def reverse_complement(s):
    ra = {65:84, 67:71, 71:67, 84:65, 78:78}
    return bytes([ra[x] for x in s[::-1]])


def hamming(a, b):
    la = len(a)
    lb = len(b)
    if la != lb:
        return la
    return sum([a[x] != b[x] for x in range(len(a))])


def get_options():
    parser = argparse.ArgumentParser(prog='process_ont1x.py')
    parser.add_argument('-i', '--fastq', help='ONT cDNA fastq')
    parser.add_argument('-p', '--prefix', help='Prefix for output files')

    options = parser.parse_args()
  
    return options

def main():
    nl = b'\n'
    dnl = b'\n+\n'
    dark = b'GGGGGGGGGGGGGGGGGGGG'
    _chunk_size = 512 # number 

    options = get_options()

    AD = b'CTACACGACGCTCTTCCGATCT'
    l = len(AD)//2
    ad1 = AD[:l]
    ad2 = AD[l:]
    l1 = len(ad1)
    l2 = len(ad2)
        
    fout1 = f'{options.prefix}_R1.fastq.gz'
    fout2 = f'{options.prefix}_R2.fastq.gz'
    foutF = f'{options.prefix}_FAIL.fastq.gz'
    raw_out1 = open(fout1, 'wb')
    raw_out2 = open(fout2, 'wb')
    raw_outF = open(foutF, 'wb')    
    fh_out1 = bgzip.BGZipWriter(raw_out1, batch_size=256)
    fh_out2 = bgzip.BGZipWriter(raw_out2, batch_size=256)
    fh_outF = bgzip.BGZipWriter(raw_outF, batch_size=256)
    r1_spool = b''
    r2_spool = b''
    F_spool = b''
    n_tot = 0
    n_fail = 0
    _spool_counter = 0

    _stream = options.fastq    
    if options.fastq == None:
        _stream = sys.stdin
    
    
    
    for entry in HTSeq.FastqReader(_stream):
        n_tot += 1
        S = entry.seq#[:60]
        Q = entry.qualstr#[:60]
        p1, p2 = S.find(ad1), S.find(ad2)
        
        if p1 == -1 and p2 == -1:
            S = reverse_complement(entry.seq)
            Q = entry.qualstr[::-1]
            p1, p2 = S.find(ad1), S.find(ad2)
            if p1 == -1 and p2 == -1:
                n_fail += 1
                nameF = bytes('@' + entry.name, encoding='ascii')
                F_spool = F_spool + nameF + nl + entry.seq + dnl + entry.qualstr + nl
                continue    
        strip_pos = -1
        if p2 - p1 == l:
            #exact match
            strip_pos = p2 + l2
        elif p1 > 0 and hamming(S[p1+l1:p1+l1+l2], ad2) < 2:
            strip_pos = p1 + l1 + l2
        elif p2 > 0 and hamming(S[p2:p2-l1], ad1) < 2:
            strip_pos = p2 + l2
        if strip_pos == -1:
            n_fail += 1
            nameF = bytes('@' + entry.name, encoding='ascii')
            F_spool = F_spool + nameF + nl + entry.seq + dnl + entry.qualstr + nl
            continue
        name1 = bytes('@' + entry.name + '/1', encoding='ascii')
        name2 = bytes('@' + entry.name + '/2', encoding='ascii')
        
        seq1 = S[strip_pos:strip_pos + 28]
        qual1 = Q[strip_pos:strip_pos + 28]
        seq2 = S[strip_pos + 28:]
        qual2 = Q[strip_pos + 28:]
        
        seq2 = reverse_complement(seq2)
        qual2 = qual2[::-1]
            
        r1_spool = r1_spool + name1 + nl + seq1 + dnl + qual1 + nl
        r2_spool = r2_spool + name2 + nl + seq2 + dnl + qual2 + nl
       
        _spool_counter += 1
        
        if _spool_counter == _chunk_size:
            fh_out1.write(r1_spool)
            fh_out2.write(r2_spool)
            fh_outF.write(F_spool)
            r1_spool = b''
            r2_spool = b''
            F_spool = b''
            _spool_counter = 0
    
    # end, write remaining spool and close files
    fh_out1.write(r1_spool)
    fh_out1.close()
    raw_out1.close()
    fh_out2.write(r2_spool)
    fh_out2.close()
    raw_out2.close()
    fh_outF.write(F_spool)
    fh_outF.close()
    raw_outF.close()

    n_pass = n_tot - n_fail
    sys.stderr.write(f"Total sequences:\t{n_tot}\n")
    f = n_fail / n_tot * 100
    sys.stderr.write(f"Failed sequences:\t{n_fail} ({f:.3f}%)\n")
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
    
        
    

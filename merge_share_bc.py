import pandas as pd
import sys
import numpy as np
import editdistance as ed
from tqdm import tqdm
import argparse
import HTSeq


# struct Barcode
# read     CCGAGCCCACGAGAC | TCGGACGATCATGGGNNNNNNNNCAAGTATGCAGCGCGCTCAAGCACGTGGATNNNNNNNNAGTCGTACGCCGATGCGAAACATCGGCCACNNNNNNNN | ATCTCGTATGCCGTCTTCTGCTTG
# template GGCTCGGGTGCTCTG | AGCCTGCTAGTACCCNNNNNNNNGTTCATACGTCGCGCGAGTTCGTGCACCTANNNNNNNNTCAGCATGCGGCTACGCTTTGTAGCCGGTGNNNNNNNN | TAGAGCATACGGCAGAAGACGAAC
# struct UMI
# AGATGTGTATAAGAGACAG | NNNNNNNNNN


def get_options():
  parser = argparse.ArgumentParser(prog='merge_share_bc.py')
  parser.add_argument('-b', '--read_bc', help='Read containing cell barcode (R2)')
  parser.add_argument('-u', '--read_umi', help='Read containing UMI (R3, RNA only)')
  parser.add_argument('-e', '--max_err', help='Max edit distance (comma separated)', nargs=3, default=[1, 1, 1])
  parser.add_argument('-d', '--debug', help='Print debug information', action='store_true')
#  parser.add_argument('-o', '--output', help='Prefix for outputfile (in pickle format)')
#  parser.add_argument('-E', '--threshold', help='General threshold for edit distance', default=1, type=int)
  
  options = parser.parse_args()
  
  return options

def main():
    options = get_options()
    
    thr = [int(x) for x in options.max_err]
    
    sp1 = 'TCGGACGATCATGGG' # [0:15]
    sp2 = 'CAAGTATGCAGCGCGCTCAAGCACGTGGAT' # [23:53]
#    sp3 = 'AGTCGTACGCCGATGCGAAACATCGGCCAC' # [61:91]
    sp3 = 'AGTCGTACGCCGATGCAAAACATAAACCAC' # [61:91]

    bc_fh = iter(HTSeq.FastqReader(options.read_bc))
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
        if e1 > thr[0] and e2 > thr[1] and e3 > thr[2]:
            continue
        comb_bc = bc_s[15:23] + bc_s[53:61] + bc_s[91:99]
        comb_ql = bc_q[15:23] + bc_q[53:61] + bc_q[91:99]
        if 'N' in comb_bc:
            continue
        
        if options.read_umi:
            umi_s = umi_r.seq.decode('ascii')
            umi_q = umi_r.qualstr.decode('ascii')
            comb_bc = comb_bc + umi_s[:10]
            comb_ql = comb_ql + umi_q[:10]
        
        sys.stdout.write(f'@{bc_r.name}\n{comb_bc}\n+\n{comb_ql}\n')
            


if __name__ == '__main__':
    main()

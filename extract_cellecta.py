import pandas as pd
import sys
import numpy as np
import editdistance as ed
from tqdm import tqdm
import argparse
import HTSeq


def get_options():
  parser = argparse.ArgumentParser(prog='extract_cellecta.py')
  parser.add_argument('-1', '--read_cellecta', help='Read containing cellecta bc (R2)')
  parser.add_argument('-2', '--read_10x', help='Read containing 10x bc (R1)')
  parser.add_argument('-w', '--whitelist', help='Whitelist from umi-tools')
  parser.add_argument('--bc14', help='Object with bc14')
  parser.add_argument('--bc30', help='Object with bc30')  
  parser.add_argument('-o', '--output', help='Prefix for outputfile (in pickle format)')
  parser.add_argument('-E', '--threshold', help='General threshold for edit distance', default=1, type=int)
  
  options = parser.parse_args()
  
  return options

def main():
    options = get_options()
    
    bc14 = pd.read_pickle(options.bc14)
    bc30 = pd.read_pickle(options.bc30)

    bc14_f = bc14.iloc[:, 2].values
    bc14_r = bc14.iloc[:, 3].values
    bc30_f = bc30.iloc[:, 2].values
    bc30_r = bc30.iloc[:, 3].values
    
    # read whitelist
    wl_a = []
    with open(options.whitelist) as fh:
        for line in fh:
            t = line.split()
            wl_a.append(t[0])

    cl_fh = iter(HTSeq.FastqReader(options.read_cellecta))
    bc_fh = iter(HTSeq.FastqReader(options.read_10x))

    cl_d = {}
    cl_df = []
    collect_cr = []
    ny = 0
    nn = 0
    for (c, b) in tqdm(zip(cl_fh, bc_fh)):
        cellecta = c.seq.decode('ascii')
        fbp1 = cellecta[:25]
        if ed.eval(fbp1, 'CCGACCACCGAACGCAACGCACGCA') > options.threshold:
            continue
        a_pos = cellecta.rfind('TGGT') 
    #    if a_pos > 40 or a_pos < 38:
    #        continue
    #    start = a_pos - 14
    #    end = a_pos + 34
        start = 25
        end = 73
        _s = cellecta[start:end]
        if ed.eval(_s[14:18],'TGGT') > 1:
            continue
    #    if len(_s) != 48:
    #        continue
        f14 = _s[:14]
        f30 = _s[-30:]
        d1 = [ed.eval(f14, x) for x in bc14_f]
        if min(d1) > options.threshold:
            continue
        d2 = [ed.eval(f30, x) for x in bc30_f]
        if min(d2) > options.threshold:
            continue
        b14 = bc14_f[np.argmin(d1)]
        b30 = bc30_f[np.argmin(d2)]
        _c = b14 + 'TGGT' + b30
        chromium = b.seq.decode('ascii')[:16]
        in_wl = 0
        bc10x = ''
        d3 = [ed.eval(chromium, x) for x in wl_a]
        if min(d3) <= options.threshold:
    #    if chromium in wl_d:
            in_wl = 1
            bc10x = wl_a[np.argmin(d3)]
        cl_df.append([_c, chromium, bc10x, in_wl])
    
    cl_df = pd.DataFrame(cl_df, columns=['cellecta', 'chromium', 'bc10x','whitelist'])
    cl_df.to_pickle(f'{options.output}.pickle')


if __name__ == '__main__':
    main()

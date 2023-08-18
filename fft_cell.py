import sys
import os
import argparse
import numpy as np
import scipy.signal as ssig
import gzip
import pandas as pd
from tqdm import tqdm


def get_options():
    parser = argparse.ArgumentParser(prog='fft_cell.py')
    parser.add_argument('-i', '--sample', help='Sample name (sample_BC_barcode.bed)', required=True)
    parser.add_argument('-p', '--path', help='Path of BED files', default='.')
    parser.add_argument('-o', '--output', help='Output file', default='output.pickle')
    parser.add_argument('-w', '--whitelist', help='Limit the analysis to these cells')
    parser.add_argument('-n', '--number', help='Limit the analysis to these much fragments', default=-1, type=int)
  
    options = parser.parse_args()
    return options


def parse_beds(sample, path, barcodes, whitelist=set(), max_lines=-1):
    csize = dict()
    for barcode in barcodes:
        filename = f'{path}/{sample}_BC_{barcode}.bed'
        if os.path.exists(filename):
            is_gzip = False
        elif os.path.exists(f'{filename}.gz'):
            is_gzip = True
            filename = f'{filename}.gz'
        else:
            sys.stderr.write(f"File type not supported for sample {sample}, {barcode}\n")
            continue    

        if is_gzip:
            fh = gzip.open(filename)
        else:
            fh = open(filename)
    
        nl = 0
        sys.stderr.write(f"Parsing {filename}")
        for line in tqdm(fh):
            nl += 1
            if nl == max_lines:
                break

            if is_gzip:
                line = line.decode('ascii')
            fields = line.split('\t')
            cell = t[3]
            frag_len = int(t[2]) - int(t[1])
        
            if not cell in csize:
                csize[cell] = []

            if len(whitelist) == 0 or cell in whitelist:
                csize[cell].append(frag_len)

    return csize                
    
def fft_cell():

    barcodes =  {'tn5':['CGTACTAG','TCCTGAGC','TCATGAGC','CCTGAGAT'],
                 'tnh':['TAAGGCGA','GCTACGCT','AGGCTCCG','CTGCGCAT']}

    d = []
    bc_tn5 = barcodes['tn5']#['CGTACTAG','TCCTGAGC','TCATGAGC','CCTGAGAT']
    bc_tnh = barcodes['tnh']#['CGTACTAG','TCCTGAGC','TCATGAGC','CCTGAGAT']

    options = get_options()

    whitelist = set()
    if options.whitelist:
        for line in open(options.whitelist):
            t = line.split()
            whitelist.add([t[0]])

    _csize_tn5 = parse_beds(options.sample, options.path, bc_tn5, whitelist, options.number)
    _csize_tnh = parse_beds(options.sample, options.path, bc_tnh, whitelist, options.number)

    csize_tn5 = dict()
    csize_tnh = dict()

    cells = set(_csize_tn5.keys()).intersection(_csize_tnh.keys())
    for cell in cells:
        csize_tn5[cell] = np.array(_csize_tn5[cell])
        csize_tnh[cell] = np.array(_csize_tnh[cell])
    
    del _csize_tn5, _csize_tnh
    
    spec_tn5 = dict.fromkeys(cells)
    spec_tnh = dict.fromkeys(cells)
    hist_tn5 = dict.fromkeys(cells)
    hist_tnh = dict.fromkeys(cells)
    
    nf = 10
    bins = np.arange(1, 2000, nf)

    for cell in tqdm(cells):
        H = np.histogram(csize_tn5[cell], bins=bins, density=True)
        a = np.fft.rfft(H[0])
        NS = len(H[0])
        xr = np.fft.rfftfreq(NS, 1/nf)
        hist_tn5[cell] = H[0]
        spec_tn5[cell] = a
        H = np.histogram(csize_tnh[cell], bins=bins, density=True)
        a = np.fft.rfft(H[0])
        NS = len(H[0])
        xr = np.fft.rfftfreq(NS, 1/nf)
        hist_tnh[cell] = H[0]
        spec_tnh[cell] = a    
        
    H_tn5 = np.array([hist_tn5[x] for x in cells])
    H_tnh = np.array([hist_tnh[x] for x in cells])
    S_tn5 = np.array([spec_tn5[x] for x in cells])
    S_tnh = np.array([spec_tnh[x] for x in cells])        
    
    df = pd.DataFrame(0, index=cells, columns=['NS_tn5', 'NS_tnH', 'NS_odds', 'NS_diff'])
    
    df['NS_tn5'] = H_tn5[:, :13].sum(1) / H_tn5[:, 13:26].sum(1)
    df['NS_tnH'] = H_tnh[:, :13].sum(1) / H_tnh[:, 13:26].sum(1)
    df['NS_odds'] = df['NS_tn5'] / df['NS_tnH']
    df['NS_diff'] = df['NS_tn5'] - df['NS_tnH']

    df.to_pickle(options.output)

if __name__ == '__main__':
    fft_cell()

import pandas as pd
import argparse
import numpy as np
import scipy.sparse as sp
import anndata as ad
from tqdm import tqdm
import gzip
import os
import sys



def get_options():
    parser = argparse.ArgumentParser(prog='parse_chromap.py')
    parser.add_argument('-i', '--sample', help='Sample name (sample_BC_barcode.bed)', required=True)
    parser.add_argument('-b', '--bin_file', help='BED file with genomic bins')
    parser.add_argument('-S', '--bin_size', help='Bin size', default=5000, type=int)
    parser.add_argument('-N', '--no_nfr', help='Discard NFR from tnH signal', action='store_true')
    parser.add_argument('-l', '--frag_len', help='Max size for NFR', default=120, type=int)
  
    options = parser.parse_args()
    return options


def parse_chromap():

    options = get_options()
    
    sample = options.sample
    bin_file = options.bin_file
    bins = pd.read_table(bin_file, header=None)
    bin_size = options.bin_size

    chroms = set(bins[0])
    offsets = dict.fromkeys(chroms)
    
    for c in chroms:
        offsets[c] = np.where(bins[0]==c)[0][0]
        
    barcodes =  {'tn5':['CGTACTAG','TCCTGAGC','TCATGAGC','CCTGAGAT'],
                 'tnH':['TAAGGCGA','GCTACGCT','AGGCTCCG','CTGCGCAT']}
                 
    whitelist = [x.split()[0] for x in open(f"{sample}_whitelist.tsv")]
    bidx = dict([(whitelist[x], x) for x in range(len(whitelist))])

    M = dict([(x, sp.lil_matrix((len(whitelist), len(bins)))) for x in barcodes.keys()])

    for tn in barcodes:
        for bc in barcodes[tn]:
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
    
            print(tn, bc)
            for line in tqdm(fh):
                if is_gzip:
                    line = line.decode('ascii')

                chrom, start, end, name, dup_count = line.split()

                if not chrom in offsets:
                    continue
                
                if tn == 'tnH' and options.no_nfr:
                    if int(end) - int(start) > options.frag_len:
                        continue
                os = int(start) // bin_size
                oe = int(end) // bin_size
                dup_count = int(dup_count)
                if os != oe:
                    c_idx = os + offsets[chrom]
                    r_idx = bidx[name]
                    M[tn][r_idx, c_idx] = M[tn][r_idx, c_idx] + 1
                    c_idx = oe + offsets[chrom]
                    r_idx = bidx[name]
                    M[tn][r_idx, c_idx] = M[tn][r_idx, c_idx] + 1
                else:
                    c_idx = os + offsets[chrom]
                    r_idx = bidx[name]
                    M[tn][r_idx, c_idx] = M[tn][r_idx, c_idx] + 1
            
    for tn in M:
        M[tn] = sp.csr_matrix(M[tn])
    
    adata = ad.AnnData(M['tn5'])
    for tn in M.keys():
        adata.layers[tn] = M[tn]
    adata.obs_names = whitelist
    adata.var_names = [f"{x[0]}:{x[1]}-{x[2]}" for x in bins.values]
    
    adata.uns['NFR'] = {'parse_NFR':options.no_nfr, 'threshold':options.frag_len}
    
    adata.write(f"{sample}.h5ad")



if __name__ == '__main__':
    parse_chromap()

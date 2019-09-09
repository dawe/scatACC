import sys
import pysam
import argparse
import scipy.sparse
import numpy as np
import pandas as pd
import pyBigWig
import scanpy as sc

def get_options():
  parser = argparse.ArgumentParser(prog='extract_bigwig.py')
  parser.add_argument('-i', '--input-file', help='Input file (.h5ad anndata)', required=True)
  parser.add_argument('-n', '--cell-number', help='Number of cell to extract (default: 0 for all)', type=int, default=0)
  parser.add_argument('-p', '--prefix', help='Prefix for output file names', default='')
  parser.add_argument('-t', '--header-tab', help='Chromosome sizes', required=True)

  options = parser.parse_args()


  return options

def get_header(f):
  header = []
  for line in open(f):
    t = line.split()
    header.append((t[0], int(t[1])))  
  return header  

def main():
  options = get_options()

  header = get_header(options.header_tab)
  chrom_sizes = dict(header)
  
  if options.prefix:
    out_prefix = "%s_" % options.prefix
  else:
    out_prefix = ''

  adata = sc.read(options.input_file)
  if options.cell_number <= 0 or options.cell_number > adata.shape[0]:
    n_extract = adata.shape[0]
  else:  
    n_extract = options.cell_number
  
  keep_cells = adata.obs.sort_values('n_regions', ascending=False).index[:n_extract]
  adata = adata[keep_cells]
  
  
  
  for cx in range(len(keep_cells)):
    fout = pyBigWig.open("%s%s.bigwig" % (out_prefix, keep_cells[cx]), 'wb')
    fout.addHeader(header)
    for _c in chrom_sizes.keys():
      cmask = np.array([x.split(':')[0] == _c for x in adata.var.index])
      if np.sum(cmask) > 0:
        fout.addEntries(_c, 0, values = adata[keep_cells[cx]][:, cmask].X, span=5000, step=5000)
    fout.close()    

if __name__ == '__main__':
  main()


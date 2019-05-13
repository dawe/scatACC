import sys
import pysam
import argparse as sp
import scipy.sparse
import numpy as np
import pandas as pd
import scanpy as sc
import pyranges as pr

def get_options():
  parser = argparse.ArgumentParser(prog='cna_analysis.py')
  parser.add_argument('-i', '--input-file', help='AnnData File with binned counts', required=True)
  parser.add_argument('-g', '--gcContent', help='GC content file (bed file)')
  parser.add_argument('-p', '--prefix', help='Prefix for output file')
  parser.add_argument('-w', '--window-size', help='Window size for CNA analysis', type=int, default=10000000)
  parser.add_argument('-s', '--step-size', help='Step size for CNA analysis', type=int, default=2000000)

  
  options = parser.parse_args()


  return options


def main():
  options = get_options()
  window_size = options.window_size
  step_size = options.step_size


  adata = sc.read(options.input_file)
  data_mat = sp.csc_matrix(adata.X) #we are going to perform lots of column slicing, CSC may be better
  gcContent = pr.read_bed(options.gcContent)
  gcContent.gcCount = gcContent.Name
  gcContent = gcContent.drop('Score').drop('Name')
  
  gcContent.data_idx = np.arange(len(gcContent))
  
  chrom_list = gcContent.Chromosome.cat.categories
  
  raw_cna = []
  raw_gc = pr.PyRanges()
  
  for _chr in chrom_list:
    chrom_size = gcContent[_chr].End.max()
    for r_start in np.arange(0, chrom_size, step_size):
      r_end = r_start + window_size
      raw_cna.append(np.array(data_mat[:, gcContent[_chr, r_start:r_end].data_idx].sum(axis=1) ).ravel())
      intv_len = window_size
      if r_end > chrom_size:
        intv_len = chrom_size - r_start
      this_gc = pr.PyRanges(chromosomes=_chr, starts=[r_start], ends=[r_end])
      this_gc.gcContent = gcContent[_chr, r_start:r_end].gcCount.sum() / intv_len
      raw_gc = pr.concat([raw_gc, this_gc])


  
  
  


if __name__ == '__main__':
  main()


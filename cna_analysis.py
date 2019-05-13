import sys
import pysam
import argparse
import scipy.sparse
import numpy as np
import pandas as pd
import scanpy as sc

def get_options():
  parser = argparse.ArgumentParser(prog='cna_analysis.py')
  parser.add_argument('-i', '--input-file', help='AnnData File with binned counts', required=True)
  parser.add_argument('-g', '--gcContent', help='GC content file (pickled dataframe)')
  parser.add_argument('-p', '--prefix', help='Prefix for output file')
  parser.add_argument('-w', '--window-size', help='Window size for CNA analysis', type=int, default=10000000)
  parser.add_argument('-s', '--step-size', help='Step size for CNA analysis', type=int, default=2000000)

  
  options = parser.parse_args()


  return options

def main():
  options = get_options()


  adata = sc.read(options.input_file)
  gcContent = pd.read_pickle(options.gcContent)
  



if __name__ == '__main__':
  main()


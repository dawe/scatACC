import sys
import argparse 
import scipy.sparse as sp
import numpy as np
import pandas as pd
import scanpy as sc
import pyranges as pr

def get_options():
  parser = argparse.ArgumentParser(prog='cna_analysis.py')
  parser.add_argument('-i', '--input-file', help='AnnData File with binned counts', nargs = '+')
  parser.add_argument('-g', '--gcContent', help='GC content file (bed file)')
  parser.add_argument('-p', '--prefix', help='Prefix for output file')
  parser.add_argument('-w', '--window-size', help='Window size for CNA analysis', type=int, default=10000000)
  parser.add_argument('-s', '--step-size', help='Step size for CNA analysis', type=int, default=2000000)
  parser.add_argument('-T', '--trim-max', help='Max copy number callable', type=int, default=6)
  parser.add_argument('--no-gc', help='Do not correct for GC content', action='store_true')

  
  options = parser.parse_args()

  if options.step_size > options.window_size:
    sys.stderr.write("Step size cannot be larger than window size\nSetting it to %d" % options.window_size)
    options.step_size = options.window_size

  return options

def gc2binidx(gcContent, resolution=0.05):
    binidx = {}
    for x in np.arange(0, 1, resolution):
        binidx[x] =  np.where((gcContent >= x  ) & (gcContent < x + resolution))[0] 
    return binidx

def main():
  options = get_options()
  window_size = options.window_size
  step_size = options.step_size


  adata = sc.read(options.input_file[0])
  data_mat = sp.csc_matrix(adata.X) #we are going to perform lots of column slicing, CSC may be better
  for f in options.input_file[1:]:
    adata = sc.read(f)
    data_mat = data_mat + sp.csc_matrix(adata.X)
  gcContent = pr.read_bed(options.gcContent)
  gcContent.gcCount = gcContent.Name
  gcContent = gcContent.drop('Score').drop('Name')
  
  gcContent.data_idx = np.arange(len(gcContent))
  
  chrom_list = gcContent.Chromosome.cat.categories
  
  raw_cna = []
  raw_gc = pr.PyRanges()
  
  nbin = 0
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
      this_gc.binidx = nbin
      raw_gc = pr.concat([raw_gc, this_gc])
      nbin += 1

  raw_cna = np.array(raw_cna) + 0.5
  binidx = gc2binidx(raw_gc.gcContent.values, resolution = 0.05)
  M_raw = np.mean(raw_cna, axis=0)

  gcR = np.arange(0, 1, 0.05)
  cna_ratio = np.zeros_like(raw_cna)
  for idx in np.arange(len(raw_cna)):
    gc_x = raw_gc.gcContent.values[idx]
    gc_bin = gcR[(gc_x >= gcR) ][-1]
    idxs = np.setdiff1d(binidx[gc_bin], [idx])
    np.random.shuffle(idxs)
    idxs = idxs[:100]
    if options.no_gc:
      cna_ratio[idx] = raw_cna[idx] / M_raw
    else:
      cna_ratio[idx] = raw_cna[idx] /  np.mean(raw_cna[idxs], axis=0)
    
  cna_calls = []
  raw_calls = []
  nbin = 0
  cna_gr = pr.PyRanges()
  for _chr in chrom_list:
    chrom_size = gcContent[_chr].End.max()
    for r_start in np.arange(0, chrom_size, window_size):
      r_end = r_start + window_size
      if r_end > chrom_size:
        r_end = chrom_size
      idxs = raw_gc[_chr, r_start:r_end].binidx.values
      cna_calls.append(np.mean(cna_ratio[idxs], axis=0))
      this_seg = pr.PyRanges(chromosomes=_chr, starts=[r_start], ends=[r_end])
      this_seg.binidx = nbin
      cna_gr = pr.concat([cna_gr, this_seg])
      nbin += 1

  idx = ["%s:%d-%d" % (x[0], x[1], x[2]) for x in cna_gr.df.values]
  cna_calls = np.array(cna_calls)      
  pd.DataFrame(cna_calls, index=idx, columns=adata.obs.index).to_pickle("%s_raw_calls.pickle" % options.prefix)
  pd.DataFrame(np.round(2 * cna_calls), index=idx, columns=adata.obs.index).to_pickle("%s_cna_calls.pickle" % options.prefix)
  


if __name__ == '__main__':
  main()


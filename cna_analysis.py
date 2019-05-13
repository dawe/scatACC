import sys
import pysam
import argparse 
import scipy.sparse as sp
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

  if options.step_size > options.window_size:
    sys.stderr.write("Step size cannot be larger than window size\nSetting it to %d" % options.window_size)
    options.step_size = options.window_size

  return options

def gc2binidx(gcContent, resolution=0.05):
  binidx = {}
  for x in np.arange(0, 1, resolution):
    binidx[x] = set(np.where((gcContent >= x  ) & (gcContent < x + resolution))[0])
  return binidx

def callcna(v, trim=4):
  v = v + 1.5
  v[v < 0] = 0
  v[v > trim] = trim
  return v // 1

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

  raw_cna = np.array(raw_cna)
  binidx = gc2binidx(raw_gc.gcContent.values, resolution = 0.05)
  cnbins = np.arange(5.5, -2.5, -1)

  
  gcR = np.arange(0, 1, 0.05)
  cna_ratio = np.zeros_like(raw_cna)
  for idx in np.arange(len(raw_cna)):
    gc_x = raw_gc.gcContent.values[idx]
    gc_bin = gcR[(idx >= gcR) ][-1]
    idxs = np.array(list(binidx[gc_bin] - {idx}))
    np.random.shuffle(idxs)
    idxs = idxs[:100]
    mean_cna = raw_cna[idxs].mean(axis=0)
    cna_r = raw_cna[idx] / mean_cna
    cna_r[np.isnan(cna_r)] = raw_cna[idx][np.isnan(cna_r)]
    cna_ratio[idx] = cna_r
    
  logcn = np.log2(cna_ratio)
  min_l = np.min(logcn[logcn > -inf])
  max_l = np.min(logcn[logcn < inf])
  logcn[logcn == -inf] = min_l
  logcn[logcn == inf] = max_l
  
  cna_calls = []
  nbin = 0
  cna_gr = pr.PyRanges()
  for _chr in chrom_list:
    chrom_size = gcContent[_chr].End.max()
    for r_start in np.arange(0, chrom_size, window_size):
      r_end = r_start + window_size
      idxs = raw_gc[_chr, r_start:r_end].binidx.values
      a_logcn = logcn[idxs].mean(axis=0)
      cna_calls.append(callcna(a_logcn, trim=4))
      this_seg = pr.PyRanges(chromosomes=_chr, starts=[r_start], ends=[r_end])
      this_seg.binidx = nbin
      cna_gr = pr.concat([cna_gr, this_seg])
      nbin += 1

  cna_calls = np.array(cna_calls)      
  
  pd.DataFrame(cna_calls, index=idx, columns=adata.obs.index).to_pickle("%s.pickle" % options.prefix)
  


if __name__ == '__main__':
  main()


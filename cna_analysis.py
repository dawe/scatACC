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
  parser.add_argument('--keep-bg', help='Keep background data (if any)', action='store_true')

  
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
  if not options.keep_bg:
    data_mat =data_mat[:-1]
      
  gcContent = pr.read_bed(options.gcContent)
  gcContent.gcCount = gcContent.Name
  gcContent = gcContent.drop('Score').drop('Name')
  
  bin_df = pd.DataFrame([x.replace(':', '\t').replace('-', '\t').split() for x in adata.var.index], columns=['Chromosome', 'Start', 'End'])
  bin_df.loc[:, 'data_idx'] = np.arange(len(bin_df))
  bin_df = pr.PyRanges(bin_df)
  
  chrom_list = gcContent.Chromosome.cat.categories
  
  raw_gc = []
  nbin = 0
  for _chr in chrom_list:
    chrom_size = gcContent[_chr].End.max()
    for r_start in np.arange(0, chrom_size, step_size):
      r_end = r_start + window_size
      intv_len = window_size
      if r_end > chrom_size:
        intv_len = chrom_size - r_start
      try:
        _gc = gcContent[_chr, r_start:r_end].gcCount.sum() / intv_len
      except IndexError:
        _gc = 0.0
      raw_gc.append([_chr, r_start, r_end, _gc, nbin])    
      nbin += 1

  raw_gc = pr.PyRanges(pd.DataFrame(raw_gc, columns = ['Chromosome', 'Start', 'End', 'gcContent', 'binidx']))
  
  raw_cna = np.zeros((len(raw_gc), data_mat.shape[0]))
  for _chr, df in raw_gc:
    for entry in df.values:
      try:
        _v = data_mat[:, bin_df[entry[0], entry[1]:entry[2]].data_idx].sum(axis=1).ravel()
      except IndexError:
        _v = 0  
      raw_cna[entry[4]] = _v

  coverage = np.array(data_mat.sum(axis=1)).ravel()
  raw_cna = np.array(raw_cna) + 0.5 # add pseudocounts
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
      cna_ratio[idx] = raw_cna[idx] / np.mean(raw_cna[idxs], axis=0)
    #  cna_ratio[idx] = raw_cna[idx] /  np.mean(raw_cna[idxs], axis=0)  
    
  cna_size = np.sum([gcContent[_chr].End.max() // window_size + 1 for _chr in chrom_list])
  cna_calls = np.zeros((cna_size, data_mat.shape[0]))
  raw_calls = np.zeros((cna_size, data_mat.shape[0]))
  nbin = 0
  cna_gr = pr.PyRanges()
  for _chr in chrom_list:
    chrom_size = gcContent[_chr].End.max()
    for r_start in np.arange(0, chrom_size, window_size):
      r_end = r_start + window_size
      if r_end > chrom_size:
        r_end = chrom_size
      idxs = raw_gc[_chr, r_start:r_end].binidx.values
      cna_calls[nbin] = np.mean(cna_ratio[idxs], axis=0)
      this_seg = pr.PyRanges(chromosomes=_chr, starts=[r_start], ends=[r_end])
      this_seg.binidx = nbin
      cna_gr = pr.concat([cna_gr, this_seg])
      nbin += 1

  idx = ["%s:%d-%d" % (x[0], x[1], x[2]) for x in cna_gr.df.sort_values('binidx').values]
  cols = adata.obs.index
  if not options.keep_bg:
    cols = cols[:-1]
  cna_calls = np.array(cna_calls)      
  pd.DataFrame(cna_calls, index=idx, columns=cols).to_pickle("%s_raw_calls.pickle" % options.prefix)
  f_corr = 2 / np.median(cna_calls) #assuming diploids
  cna_calls = np.round(f_corr * cna_calls)
  cna_calls[cna_calls > options.trim_max] = options.trim_max
  pd.DataFrame(cna_calls, index=idx, columns=cols).to_pickle("%s_cna_calls.pickle" % options.prefix)
  


if __name__ == '__main__':
  main()


import pyBigWig as pbw
import numpy as np
import sys
import argparse
import os
import pywt

def get_options():
  parser = argparse.ArgumentParser(prog='block_deconvolve.py')
  parser.add_argument('-p', '--prefix', help='Prefix for all files', default='out')
  parser.add_argument('-s', '--stepsize', help='Bin size for output', default=1000, type=int)
  parser.add_argument('-5', '--tn5', help='BigWig file with tn5 signal')
  parser.add_argument('-H', '--tnH', help='BigWig file with tnH signal')
  parser.add_argument('-D', '--write-diff', help='Write difference file', action='store_true')
  parser.add_argument('-S', '--scale', help='Scaling factor after normalization', default=1e4, type=float)
  parser.add_argument('-A', '--asinh', help='Apply asinh transformation to diff file', action='store_true')
  parser.add_argument('-M', '--smooth', help='Apply wavelet smoothing up to this level', default=0)  
  
  options = parser.parse_args()
  
  return options

def wavelet_smooth(data, level, wavelet='sym7'):
  coefs = pywt.wavedec(data, wavelet)
  level = int(level)
  if level >= len(coefs):
    level = len(coefs) - 1
  if level == 0:
    return data
  for l in range(level):
    nl = l + 1
    coefs[-nl] = np.zeros_like(coefs[-nl])
  return pywt.waverec(coefs, wavelet)  

def deconvolve():
  options = get_options()
  # tn5
  ftn5 = f'{options.prefix}_tn5.bigwig'
  if not os.path.exists(ftn5):
    ftn5 = options.tn5
  tn5 = pbw.open(ftn5, 'rb')
  
  ftnH = f'{options.prefix}_tnH.bigwig'
  if not os.path.exists(ftnH):
    ftnH = options.tnH
  tnH = pbw.open(ftnH, 'rb')

  stepsize=options.stepsize
  chroms = tn5.chroms()
  if options.write_diff:
    rat_diff = pbw.open(f'{options.prefix}_diff.bigwig', 'wb')
    rat_diff.addHeader([(x, chroms[x]) for x in chroms])
  rat_5 = pbw.open(f'{options.prefix}_dcn_tn5.bigwig', 'wb')
  rat_5.addHeader([(x, chroms[x]) for x in chroms])
  rat_H = pbw.open(f'{options.prefix}_dcn_tnH.bigwig', 'wb')
  rat_H.addHeader([(x, chroms[x]) for x in chroms])

  for chrom in chroms:
    print(chrom)
    n_bins = chroms[chrom] // stepsize
    v1 = np.array(tn5.stats(chrom, nBins=n_bins))
    v2 = np.array(tnH.stats(chrom, nBins=n_bins))
    v1[np.isnan(v1)] = 0
    v2[np.isnan(v2)] = 0
    v1 = wavelet_smooth(v1, options.smooth)
    v2 = wavelet_smooth(v2, options.smooth)
    v1 = v1 / np.linalg.norm(v1) * options.scale
    v2 = v2 / np.linalg.norm(v2) * options.scale
    V = v1 - v2
    V[np.isnan(V)] = 0
    if options.write_diff:
      if options.asinh:
        V = np.arcsinh(V)
      rat_diff.addEntries(chrom, 1, values=V, span=stepsize, step=stepsize)
    Vd=np.zeros_like(V)
    Vd[V > 0] = V[V > 0]
    rat_5.addEntries(chrom, 1, values=Vd, span=stepsize, step=stepsize)
    Vd=np.zeros_like(V)
    Vd[V < 0] = V[V < 0]
    rat_H.addEntries(chrom, 1, values=-Vd, span=stepsize, step=stepsize)
  if options.write_diff:
    rat_diff.close()
  rat_5.close()
  rat_H.close()

if __name__ == '__main__':
  deconvolve()

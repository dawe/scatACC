import sys
import pysam
import argparse
import scipy.sparse
import numpy as np
import anndata
from joblib import delayed, Parallel
import pandas as pd


def get_options():
  parser = argparse.ArgumentParser(prog='peak_counter.py')
  parser.add_argument('-p', '--peaks-file', help='Any BED-like file with peaks', required=True)
  parser.add_argument('-b', '--bamfile', help='BAM file with alignments', required=True)
  parser.add_argument('-o', '--output', help='Prefix for output file', default='output')
  parser.add_argument('-q', '--quality', help='Min. quality for an alignment to be parsed', type=int, default=15)
  parser.add_argument('-B', '--binary', help='Ouput binary counts over intervals instead of counts', action='store_true')
  parser.add_argument('-t', '--threads', help='Number of processing threads', type=int, default=1)
  
  
  options = parser.parse_args()

  return options

def list_cells(header):
  return list(set([x['SM'] for x in header['RG']]))
  
def get_regions(intervals):
  return [f'{x[0]}:{x[1]}-{x[2]}' for x in intervals]

def get_intervals(peaks_file):
  intvl_tuples = [x.split()[:3] for x in open(peaks_file) if not x.startswith('#')]
  return [(x[0], int(x[1]), int(x[2])) for x in intvl_tuples]


def count_intervals(chunk_intvl, bc_list, bam_file, n_regions, quality):
  bam_in = pysam.Samfile(bam_file, 'rb')
  n_cells = len(bc_list)
  counter = dict([(x, 0) for x in bc_list])
  chunk_matrix = scipy.sparse.lil_matrix((n_regions, n_cells), dtype=int)
  id2sm = dict([(x['ID'], x['SM']) for x in bam_in.header['RG']])
  for chrom, start, end, nl in chunk_intvl:
#    sys.stderr.write(f'working on {chrom}:{start}-{end}, {nl} ')
    if not chrom in bam_in.references:
      continue
    for alignment in bam_in.fetch(chrom, start, end):
      if alignment.is_proper_pair and alignment.mapq >= quality and alignment.is_read1 and not alignment.is_duplicate:
        counter[id2sm[alignment.get_tag('RG')]] += 1
    chunk_matrix[nl] = [counter[x] for x in bc_list]
    counter = dict([(x,0) for x in bc_list])
#    sys.stderr.write(" done\n")
  bam_in.close()
  return chunk_matrix


def main():
  options = get_options()

  prefix = options.output
  # we assume all barcodes are stored as samples in the input bamfile header (by RG)
  bam_in = pysam.Samfile(options.bamfile, 'rb')
  bc_list = list_cells(bam_in.header)
  bam_in.close()
  
#  counter = dict([(x, 0) for x in bc_list])

  intervals = get_intervals(options.peaks_file)
  regions = get_regions(intervals)
  N_regions = len(regions)
  N_cells = len(bc_list)
#  dtype = np.uint32 #np.uint16 shoudl be enough, but you neve know
  
#  if options.binary:
#    dtype = np.bool  
  count_matrix = scipy.sparse.lil_matrix((N_regions, N_cells), dtype=int)
    
#  id2sm = dict([(x['ID'], x['SM']) for x in bam_in.header['RG']])
  
  intervals = [(*intervals[x], x) for x in range(len(intervals))]
  chunk_l = len(intervals) // options.threads
  chunks = [intervals[(chunk_l*x):(chunk_l*x+chunk_l)] for x in range(options.threads)]
  
#  count_matrix = count_intervals(intervals, bc_list, options.bamfile, N_regions, options.quality)
  
  chunk_counts = Parallel(n_jobs=options.threads)(delayed(count_intervals)(i, bc_list, options.bamfile, N_regions, options.quality) for i in chunks)
  count_matrix = np.sum(chunk_counts)
#  for chrom, start, end, nl in intervals:
#    if not chrom in bam_in.header.references:
#      continue
#    for alignment in bam_in.fetch(chrom, start, end):
#      if alignment.is_proper_pair and alignment.mapq >= options.quality and alignment.is_read1:
#        if options.binary:
#          counter[id2sm[alignment.get_tag('RG')]] = 1
#        elif not alignment.is_duplicate:
#          counter[id2sm[alignment.get_tag('RG')]] += 1
#      count_matrix[nl] = [counter[x] for x in bc_list]
#    counter = dict([(x,0) for x in bc_list])

#  bam_in.close()
    
  fout = f'{prefix}.h5ad'
  count_matrix = scipy.sparse.csr_matrix(count_matrix.T) #convert
  n_cells = pd.DataFrame(np.array(np.sum(count_matrix > 0, axis=0)).ravel(), index=regions, columns=['n_cells'])
  n_regions = pd.DataFrame(np.array(np.sum(count_matrix > 0, axis=1)).ravel(), index=bc_list, columns=['n_regions'])
  adata = anndata.AnnData(count_matrix, obs=n_regions, var=n_cells)
  if options.binary:
    adata.X = (adata.X > 0).astype(int)
  adata.write_h5ad(fout)
    

if __name__ == '__main__':
  main()


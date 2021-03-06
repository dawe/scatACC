import sys
import pysam
import argparse
import scipy.sparse
import numpy as np
import pandas as pd

def get_options():
  parser = argparse.ArgumentParser(prog='prom_act.py')
  parser.add_argument('-g', '--genes-file', help='6 columns BED file with gene name in 4th and strand in 6th', required=True)
  parser.add_argument('-b', '--bamfile', help='BAM file with alignments', required=True)
  parser.add_argument('-o', '--output', help='Prefix for output file')
  parser.add_argument('-q', '--quality', help='Min. quality for an alignment to be parsed', type=int, default=15)
  parser.add_argument('-D', '--decay', help='Size of decay factor', default=1000)
  parser.add_argument('-w', '--window', help='Size of the window around the promoter', default=5000)

  
  options = parser.parse_args()


  return options

def list_cells(header):
  return list(set([x['SM'] for x in header['RG']]))
  
def get_regions(peaks_file):
  region_tuples = [x.split()[:3] for x in open(peaks_file) if not x.startswith('#')]
  return ["%s:%s-%s" % (x[0], x[1], x[2]) for x in region_tuples]

def main():
  options = get_options()
  
  win_size = options.window
  decay = options.decay


  # we assume all barcodes are stored as samples in the input bamfile header (by RG)
  bam_in = pysam.Samfile(options.bamfile, 'rb')
  bc_list = list_cells(bam_in.header)
  
  regions = pd.read_table(options.genes_file, header=None, comment='#')
  
  N_cells = len(bc_list)
  N_regions = len(regions)
  
  act_matrix = scipy.sparse.lil_matrix((N_regions, N_cells), dtype=float)
  
  id2sm = dict([(x['ID'], x['SM']) for x in bam_in.header['RG']])
  id2idx = dict([(bc_list[x], x) for x in range(len(bc_list))])


  D = np.abs(np.arange(2*win_size) - win_size) / decay #1000 should be a parameter
  Ed = np.exp(-D)
  
  region_names = []
  for idx in regions.index:
    chrom, start, end, g_name, _, g_strand = regions.loc[idx]
    if g_strand == '-':
      start, end = end, start
    f_start = start - win_size
    f_end = start + win_size
    if f_start < 0:
      f_start = 0
    if f_end > bam_in.header.get_reference_length(chrom):
      f_end = bam_in.header.get_reference_length(chrom)
    region_names.append("%s:%d-%d:%s" % (chrom, start, end, g_name))
    count_matrix = np.zeros((N_cells, 2*win_size))
    for alignment in bam_in.fetch(chrom, f_start, f_end):
      if alignment.is_proper_pair and alignment.mapq >= options.quality and alignment.is_read1 and not alignment.is_duplicate:
        offset = alignment.pos - start + win_size
        cell_n = id2idx[id2sm[alignment.get_tag('RG')]]
        count_matrix[cell_n, offset] += 1
    act_matrix[idx] = np.sum(count_matrix * Ed, axis=1)


  fout = "%s.npz" % options.output
  act_matrix = scipy.sparse.csr_matrix(act_matrix.T) #convert
  np.savez(fout, data = act_matrix.data, indices=act_matrix.indices, 
             indptr=act_matrix.indptr, shape=act_matrix.shape, bc_list=bc_list, 
             region_names=np.array(region_names))


if __name__ == '__main__':
  main()


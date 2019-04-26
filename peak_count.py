import sys
import pysam
import argparse
import scipy.sparse
import numpy as np
import anndata

def get_options():
  parser = argparse.ArgumentParser(prog='peak_counter.py')
  parser.add_argument('-p', '--peaks-file', help='Any BED-like file with peaks', required=True)
  parser.add_argument('-b', '--bamfile', help='BAM file with alignments', required=True)
  parser.add_argument('-o', '--output', help='Prefix for output file')
  parser.add_argument('-q', '--quality', help='Min. quality for an alignment to be parsed', type=int, default=15)
  parser.add_argument('-B', '--output-bed', help='Output matrix as bed file (with header)', action='store_true')
  parser.add_argument('-S', '--output-sparse', help='Output matrix as sparse matrix', action='store_true')
  parser.add_argument('-A', '--output-anndata', help='Output matrix as AnnData, suitable for scanpy', action='store_true')
  parser.add_argument('-c', '--output-counts', help='Ouput counts over intervals instead of binary data', action='store_true')
  
  options = parser.parse_args()
  if not options.output_bed and not options.output_sparse and not options.output_anndata:
    # for the time being, default output is bed
    options.output_bed = True

  return options

def list_cells(header):
  return [x['SM'] for x in header['RG']]
  
def get_regions(peaks_file):
  region_tuples = [x.split()[:3] for x in open(peaks_file) if not x.startswith('#')]
  return ["%s:%s-%s" % (x[0], x[1], x[2]) for x in region_tuples]

def main():
  options = get_options()

  if options.output:
    prefix = options.output
    if options.output_bed:
      fout = open("%s.bed" % (prefix), 'w')
  else:
    fout = sys.stdout
    prefix = 'output'
    
  # we assume all barcodes are stored as samples in the input bamfile header (by RG)
  bam_in = pysam.Samfile(options.bamfile, 'rb')
  bc_list = list_cells(bam_in.header)
  
  counter = dict([(x, 0) for x in bc_list])

  if options.output_sparse or options.output_anndata:
    regions = get_regions(options.peaks_file)
    N_regions = len(regions)
    N_cells = len(bc_list)
    if options.output_counts:
      dtype = np.uint32 #np.uint16 shoudl be enough, but you neve know
    else:
      dtype = np.bool  
    count_matrix = scipy.sparse.lil_matrix((N_regions, N_cells), dtype=dtype)
    
  id2sm = dict([(x['ID'], x['SM']) for x in bam_in.header['RG']])
  
  if options.output_bed:
    spool = ['chrom', 'start', 'end'] + bc_list
    fout.write("\t".join(spool) + "\n")

  for nl, line in enumerate(open(options.peaks_file)):
    # TODO: use some better interval reader
    if line.startswith('#'):
      continue
    chrom, start, end = line.split()[:3]
    for alignment in bam_in.fetch(chrom, int(start), int(end)):
      if alignment.is_proper_pair and alignment.mapq >= options.quality and alignment.is_read1:
        if options.output_counts:
          # TODO: deduplicate on the fly by hashing reads
          counter[id2sm[alignment.get_tag('RG')]] += 1
        else:  
          counter[id2sm[alignment.get_tag('RG')]] = 1
    if options.output_bed:
      spool = [chrom, start, end] + ["%d" % counter[x] for x in bc_list]
      fout.write("\t".join(spool) + "\n")
    if options.output_sparse or options.output_anndata:
      count_matrix[nl] = [counter[x] for x in bc_list]
    counter = dict([(x,0) for x in bc_list])

  bam_in.close()
  if options.output_sparse:
    fout = "%s.npz" % prefix
    count_matrix = scipy.sparse.csr_matrix(count_matrix) #convert
    np.savez(fout, data = count_matrix.data, indices=count_matrix.indices, 
             indptr=count_matrix.indptr, shape=count_matrix.shape, bc_list=bc_list)

    # this to be load as following:
    # loader = np.load(file.npz)
    # count_matrix = scipy.sparse.csr_matrix((loader['data'], loader['indices'], loader['indptr']), shape=loader['shape'])             
  if options.output_anndata:
    import pandas as pd
    fout = "%s.h5ad" % prefix
    count_matrix = scipy.sparse.csr_matrix(count_matrix.T) #convert
    n_cells = pd.DataFrame(np.array(np.sum(count_matrix > 0, axis=0)).ravel(), index=regions, columns=['n_cells'])
    n_regions = pd.DataFrame(np.array(np.sum(count_matrix > 0, axis=1)).ravel(), index=bc_list, columns=['n_regions'])
    adata = anndata.AnnData(count_matrix, obs=n_regions, var=n_cells)
    adata.write_h5ad(fout)
    

if __name__ == '__main__':
  main()


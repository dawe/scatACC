import sys
import pysam
import argparse
import scipy.sparse

def get_options():
  parser = argparse.ArgumentParser(prog='peak_counter.py')
  parser.add_argument('-p', '--peaks-file', help='Any BED-like file with peaks', required=True)
  parser.add_argument('-b', '--bamfile', help='BAM file with alignments', required=True)
  parser.add_argument('-o', '--output', help='Prefix for output file')
  parser.add_argument('-q', '--quality', help='Min. quality for an alignment to be parsed', type=int, default=15)
  parser.add_argument('-B', '--output-bed', help='Output matrix as bed file (with header)', action='store_true')
  parser.add_argument('-S', '--output-sparse', help='Output matrix as sparse matrix', action='store_true')
  parser.add_argument('-A', '--output-ann', help='Output matrix as AnnData, suitable for scanpy', action='store_true')
  parser.add_argument('-c', '--output-counts', help='Ouput counts over intervals instead of binary data', action='store_true')
  
  options = parser.parse_args()
  if not options.output_bed and not options.output_sparse and not options.output_ann:
    # for the time being, default output is bed
    options.output_bed = True

  return options

def list_cells(header):
  return [x['SM'] for x in header['RG']]
  

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

  if options.output_sparse:
    N_regions = len([x for x in open(options.peaks_file) if not x.startswith('#')])
    N_cells = len(bc_list)
    count_matrix = scipy.sparse.lil_matrix((N_regions, N_cells), dtype=int)
    
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
          counter[id2sm[dict(alignment.tags)['RG']]] += 1
        else:  
          counter[id2sm[dict(alignment.tags)['RG']]] = 1
    if options.output_bed:
      spool = [chrom, start, end] + ["%d" % counter[x] for x in bc_list]
      fout.write("\t".join(spool) + "\n")
      counter = dict([(x,0) for x in bc_list])
    if options.output_sparse:
      count_matrix[nl] = [counter[x] for x in bc_list]

  bam_in.close()
  if output_sparse:
    fout = "%s.npz" % prefix
    count_matrix = scipy.sparse.csc_matrix(count_matrix) #covert
    np.savez(fout, data = count_matrix.data, indices=count_matrix.indices, 
             indptr=count_matrix.indptr, shape=count_matrix.shape, bc_list=bc_list)

    # this to be load as following:
    # loader = np.load(file.npz)
    # count_matrix = scipy.sparse.csc_matrix((loader['data'], loader['indices'], loader['indptr']), shape=loader['shape'])             

if __name__ == '__main__':
  main()


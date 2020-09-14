import pysam
import pandas as pd
import sys
import argparse
import os

def get_options():
  parser = argparse.ArgumentParser(prog='bc2rg.py')
  parser.add_argument('-c', '--clusters', help='File with cluster assignments (either tsv or pickle)', required=True)
  parser.add_argument('-b', '--bamfile', help='BAM file with alignments', required=True)
  parser.add_argument('-o', '--output', help='Prefix for output BAM files')
  
  options = parser.parse_args()
  
  return options

def main():
  options = get_options()
  
  bam_in = pysam.Samfile(options.bamfile, 'rb')
  
  header_in = bam_in.header.as_dict()
  
  if not options.output:
    prefix, _ = os.path.splitext(os.path.basename(options.bamfile))
  else:
    prefix = 'output'

  try:
    clusters = pd.read_table(options.clusters, index_col=0, header=None)
  except:
    clusters = pd.read_pickle(options.clusters)

  merged = clusters.astype(str).agg('_'.join, axis=1)     
  
  cell2group = dict(merged)
  groups = set(cell2group.values())
  
  header_out = header_in.copy()
  header_out['RG'] = []
  headers = dict([(g, header_out.copy()) for g in groups])
  bamnames = dict([(g, f'{prefix}_{g}.bam') for g in groups])
  
  cell_to_remove = set()
  for rid in header_in['RG']:
    cell_name =  rid['ID']
    if not cell_name in cell2group:
        cell_to_remove.add(cell_name)
    else:      
      g = str(cell2group[cell_name])
      rid['SM'] = g
      headers[g]['RG'].append(rid)

  bam_out = dict([(g, pysam.Samfile(bamnames[g], 'wb', header=headers[g])) for g in groups])

  
  for r in bam_in:
    rg = dict(r.tags)['RG'] 
    if rg in cell_to_remove:
      continue
    g = cell2group[rg]
    bam_out[g].write(r)

  bam_in.close()
  [f.close() for f in bam_out]    

if __name__ == '__main__':
  main()
  

import pysam
import pandas as pd
import sys
import argparse

def get_options():
  parser = argparse.ArgumentParser(prog='bc2rg.py')
  parser.add_argument('-c', '--clusters', help='File with cluster assignments (either tsv or pickle)', required=True)
  parser.add_argument('-b', '--bamfile', help='BAM file with alignments', required=True)
  parser.add_argument('-o', '--output', help='Output BAM file', default='output.bam')
  parser.add_argument('-O', '--stdout', help='Output to stdout', action='store_true')
  parser.add_argument('-I', '--stdin', help='Input from stdin', action='store_true')
  parser.add_argument('-B', '--send-bg', help='Send unmatched cells into background (default remove)', action='store_true')
  
  options = parser.parse_args()
  
  return options

def main():
  options = get_options()
  
  if options.stdin:
    bam_in = pysam.Samfile('-','r')
  else:    
    bam_in = pysam.Samfile(options.bamfile, 'rb')
  
  header_in = bam_in.header.as_dict()

  try:
    clusters = pd.read_table(options.clusters, index_col=0, header=None)
  except UnicodeDecodeError:
    clusters = pd.read_pickle(options.clusters)

  clusters = dict(clusters.iloc[:, 0])
  header_out = header_in.copy()
  header_out['RG'] = []
  cell_to_remove = set()
  for rid in header_in['RG']:
    cell_name =  rid['ID']
    if not cell_name in clusters:
      if options.send_bg:
        rid['SM'] = 'Background'
        header_out['RG'].append(rid)
      else:
        cell_to_remove.add(cell_name)
    else:      
      rid['SM'] = str(clusters[cell_name])
      header_out['RG'].append(rid)

  if options.stdout:
    bam_out = pysam.Samfile("-", 'wb', header=header_out)
  else:
    bam_out = pysam.Samfile(options.output, 'wb', header=header_out)
  
  for r in bam_in:
    rg = dict(r.tags)['RG'] 
    if rg in cell_to_remove:
      continue
    bam_out.write(r)

  bam_in.close()
  bam_out.close()    

if __name__ == '__main__':
  main()
  

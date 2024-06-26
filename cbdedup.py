import pysam
import argparse
import sys
from collections import Counter
from natsort import natsorted
import os

def get_options():
  parser = argparse.ArgumentParser(prog='cbdedup.py')
  parser.add_argument('-i', '--bamfile', help='BAM file with alignments')
  parser.add_argument('-o', '--output', help='Output BAM file', default='output.bam')
  parser.add_argument('-B', '--buffer', help='Buffer size for dedup', default=65536)
  parser.add_argument('-O', '--stdout', help='Output to stdout', action='store_true')
  parser.add_argument('-I', '--stdin', help='Input from stdin', action='store_true')
  
  options = parser.parse_args()
  
  return options


def hash_read(r):
  s1 = r.reference_start
  c1 = r.reference_id
  rg = r.get_tag('RG')
  cigar = r.cigarstring
  if r.is_paired:
    s2 = r.next_reference_start
    c2 = r.next_reference_id
  else:
    s2 = -1
    c2 = -1
  hash = "%d:%d:%d:%d:%s:%s" % (s1, c1, s2, c2, rg, cigar)
  return hash

def main():
  
  options = get_options()
  
  if options.stdin:
    bam_in = pysam.AlignmentFile('-', 'rb')
  else:
    bam_in = pysam.AlignmentFile(options.bamfile, 'rb')

  if options.stdout:
    stat_file = 'stdout_dup_stats.txt'
    options.output = sys.stdout
  else:
    stat_file = options.output.replace('.bam', '') + '_dup_stats.txt'

  bam_out = pysam.AlignmentFile(options.output, 'wb', template=bam_in)

  buffer = []
  hash_in_buffer = set()
  hash_counter = Counter()

  for alignment in bam_in:
    alignment.is_duplicate = False #undedup if any
    h_ali = hash_read(alignment)
    if h_ali in hash_in_buffer and not alignment.is_unmapped:
          alignment.is_duplicate = True
          hash_counter.update([h_ali]) 
    hash_in_buffer.add(h_ali)
    buffer.append(alignment)
    if len(buffer) == options.buffer:
      for ali_out in buffer:
        bam_out.write(ali_out)
      buffer = []
      hash_in_buffer = set()

  #clean buffer at the end
  for ali_out in buffer:
    bam_out.write(ali_out)

  dup_count = Counter(hash_counter.values())
  with open(stat_file, 'w') as fout:
    fout.write('clone_size\tn_reads\n')
    for k in natsorted(dup_count.keys()):
      fout.write(f"{k}\t{dup_count[k]}\n")

  bam_out.close()
  bam_in.close()


if __name__ == '__main__':
  main()

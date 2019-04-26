import pysam
import argparse

BUFLEN = 8192

def get_options():
  parser = argparse.ArgumentParser(prog='cbdedup.py')
  parser.add_argument('-i', '--bamfile', help='BAM file with alignments', required=True)
  parser.add_argument('-o', '--output', help='Output BAM file', default='output.bam')
  
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
  
  bam_in = pysam.Samfile(options.bamfile, 'rb')
  bam_out = pysam.Samfile(options.output, 'wb', template=bam_in)

  buffer = []
  hash_in_buffer = set()

  for alignment in bam_in:
    duplicate = False
    h_ali = hash_read(alignment)
    if 	h_ali in hash_in_buffer:
  	  alignment.is_duplicate = True
    hash_in_buffer.add(h_ali)
    buffer.append(alignment)
    if len(buffer) == BUFLEN:
      for ali_out in buffer:
        bam_out.write(ali_out)
      buffer = []
      hash_in_buffer = set()	

  #clean buffer at the end
  for ali_out in buffer:
    bam_out.write(ali_out)


  bam_out.close()
  bam_in.close()


if __name__ == '__main__':
  main()
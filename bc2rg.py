import HTSeq
import pysam
import sys
import argparse

def get_options():
  parser = argparse.ArgumentParser(prog='bc2rg.py')
  parser.add_argument('-w', '--whitelist', help='Whitelist file from UMI_tools', required=True)
  parser.add_argument('-i', '--barcodes', help='Fastq file with barcodes', required=True)
  parser.add_argument('-b', '--bamfile', help='BAM file with alignments', required=True)
  parser.add_argument('-o', '--output', help='Output BAM file', default='output.bam')
  parser.add_argument('-G', '--group', help='Additional tag for read group definition', default = '')
  parser.add_argument('-k', '--keep-unmatched', help='Keep barcodes not in whitelist as background', action='store_true')
  parser.add_argument('-O', '--stdout', help='Output to stdout', action='store_true')
  parser.add_argument('-I', '--stdin', help='Input from stdin', action='store_true')
  
  options = parser.parse_args()
  
  return options

def build_coder(wl_file):
  coder = {}
  for line in open(wl_file):
    t = line.split()
    coder[t[0]] = t[0]
    if len(t) == 3:
      for bc in t[1].split(','):
        coder[bc] = t[0]
        
  return coder

def read_barcodes(bc_file, coder):
  namer = {}
  for r in HTSeq.FastqReader(bc_file):
    if r.seq.decode() in coder:
      namer[r.name.split()[0]] = coder[r.seq.decode()]
  
  return namer    

def build_header(in_header, coder, group = '', keep_unmatched=True):
  header = dict()
  # default keys
  for k in ['HD', 'SQ', 'PG']:
    if k in in_header:
      header[k] = in_header[k]
  
  # read groups
  if 'RG' in in_header:
    _pl = in_header['RG'][0]['PL']
    _pu = in_header['RG'][0]['PU']
    _lb = in_header['RG'][0]['LB']
    _cn = in_header['RG'][0]['CN']
  else:
    _pl = _pu = _lb = _cn = 'default'
  _rg = []
  for bc in coder.values():
    _rg.append({'ID': "%s_%s" % (bc, group),
                'PL':_pl,
                'PU':_pu,
                'LB':_lb,
                'SM': bc,
                'CN':_cn })    
  if keep_unmatched:
    _rg.append({'ID': "Background_%s" % group,
                'PL':_pl,
                'PU':_pu,
                'LB':_lb,
                'SM': "Background",
                'CN':_cn })    
  header['RG'] = _rg
  return header                
  

def main():
  options = get_options()
  
  sys.stderr.write("Reading Whitelist\n")  
  coder = build_coder(options.whitelist)

  sys.stderr.write("Reading Fastq file with barcodes\n")  
  namer = read_barcodes(options.barcodes, coder)

  sys.stderr.write("Processing BAM File\n")  
  
  if options.stdin:
    bam_in = pysam.Samfile('-','r')
  else:    
    bam_in = pysam.Samfile(options.bamfile, 'rb')
  
  header = build_header(bam_in.header, coder, options.group, options.keep_unmatched)
  
  if options.stdout:
    bam_out = pysam.Samfile("-", 'wb', header=header)
  else:
    bam_out = pysam.Samfile(options.output, 'wb', header=header)
  
  for r in bam_in:
    if r.qname in namer:
      rg = "%s_%s" % (namer[r.qname], options.group)
    else:
      if not options.keep_unmatched:
        continue
      rg = 'Background_%s' % options.group

    tags = [x for x in r.tags if x[0] != 'RG']
    tags.append(('RG', rg))
    r.tags = tags
    bam_out.write(r)

  bam_in.close()
  bam_out.close()    

if __name__ == '__main__':
  main()

import sys
import pysam
import argparse

def get_options():
  parser = argparse.ArgumentParser(prog='CBfromRG.py')
  parser.add_argument('-i', '--bam_in', help='Input BAM file', required=True)
  parser.add_argument('-o', '--bam_out', help='Output BAM file', required=True)
  parser.add_argument('-b', '--background', help='Include background', action='store_true')


def main():
   options = get_options()
   
    bam_in = pysam.AlignmentFile(options.bam_in, 'rb')
    in_header = bam_in.header
    
    bam_out = pysam.AlignmentFile(options.bam_out, 'wb', header=in_header)
    rg_sm = {}
    
    for record in in_header['RG']:
        rg_sm[record['ID']] = record['SM']
    
    for record in bam_in:
        tags = record.tags.copy()
        rg = [x[1] for x in tags if x[0]=='RG'][0]
        sm = rg_sm[rg]
        if not options.background and sm != 'Background':
    	    tags.append(('CB', sm))
    	    record.tags = tags
        bam_out.write(record)
    
    bam_in.close()
    bam_out.close()    

if __name__ == '__main__':
  main()

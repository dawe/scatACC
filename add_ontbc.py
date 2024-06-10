import pysam
import sys
from tqdm import tqdm
import argparse
#from umi_tools import UMIClusterer


def hamming(a, b):
    la = len(a)
    lb = len(b)
    if la != lb:
        return la
    return sum([a[x] != b[x] for x in range(len(a))])

def find_bc(CB, L, thr=1):
    pos = [x for x in range(len(L)) if hamming(L[x], CB) <= thr]
    if len(pos) > 0:
        return (L[pos[0]], 1)
    return (CB, 0)
    

def get_cb_struct(S):
    c, u = S.split(',')
    _, cb_start, cb_len = c.split(':')
    _, umi_start, umi_len = u.split(':')    
    cb_start = int(cb_start)
    cb_end = cb_start + int(cb_len)
    umi_start = int(umi_start)
    umi_end = umi_start + int(umi_len)
    return cb_start, cb_end, umi_start, umi_end

def get_options():
    parser = argparse.ArgumentParser(prog='process_ont1x.py')
    parser.add_argument('-f', '--fastq', help='Fastq file containing barcodes', required=True)
    parser.add_argument('-b', '--bamfile', help='BAM file containing alignments')
    parser.add_argument('-w', '--whitelist', help='Whitelist file')
    parser.add_argument('-o', '--output', help='Output BAM file')
    parser.add_argument('-S', '--struct', help='Barcode Structure', default='C:0:16,U:16:12')

    options = parser.parse_args()
  
    return options


def main():

    options = get_options()
    correct_bc = False
    cb_start, cb_end, umi_start, umi_end = get_cb_struct(options.struct)    
    
    
    wl_list = []
    
    if options.whitelist:
        for line in open(options.whitelist):
            t = line.split()
            wl_list.append(t[0])
            
    if len(wl_list) > 0:
        correct_bc = True

    bamfile = options.bamfile
    outfile = options.output
    mode_in = 'rb'
    mode_out = 'wb'
    
    
    if not options.bamfile:
        bamfile = '-'
        mode_in = 'r'
    if not options.output:
        outfile = '-'
        mode_out = 'w'

    fq_reader = pysam.FastqFile(options.fastq)
    bam_reader = pysam.AlignmentFile(bamfile, mode_in)
    sam_writer = pysam.AlignmentFile(outfile, mode_out, template=bam_reader)
    
    I = iter(zip(fq_reader, bam_reader))
    
    for seq, ali in I:
        r1seq = seq.sequence
        cb = r1seq[cb_start:cb_end]
        umi = r1seq[umi_start:umi_end]
        
        ali.tags = ali.tags + [('CR', cb),('UR',umi)]
                
        if correct_bc:
            c_cb, found = find_bc(cb, wl_list)
            if found:
                ali.tags = ali.tags + [('CB', c_cb),('UB',umi)]
        
        sam_writer.write(ali)

if __name__ == '__main__':
    main()

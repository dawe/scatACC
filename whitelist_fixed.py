import sys
import argparse
import pysam
from collections import Counter

def get_options():
    parser = argparse.ArgumentParser(prog='whitelist_fixed.py')
    parser.add_argument('-i', '--input_reads', help='File containing cell barcodes')
    parser.add_argument('-c', '--cells', help='Number of cells', default=5000, type=int)
    parser.add_argument('-r', '--reads', help='Number of reads', default=100000000, type=int)

    options = parser.parse_args()
  
    return options


whitelist = Counter()

def parse_whitelist():
    options = get_options()
    
    _MAXREADS = options.reads
    _NCELLS = options.cells
    
    nl = 0
    I = iter(pysam.FastqFile(options.input_reads, 'rb'))
    while nl <= _MAXREADS:
        entry = next(I)
        whitelist.update([entry.sequence])
        nl += 1
    
    for bc in whitelist.most_common(_NCELLS):
        sys.stdout.write(f'{bc[0]}\n')
        
    




if __name__ == '__main__':
    parse_whitelist()

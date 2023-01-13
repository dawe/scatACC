import sys
import yaml
import numpy as np
from bx.intervals import (Intersecter, Interval)
import anndata as ad
import pandas as pd
from tqdm import tqdm
import argparse
import scipy.sparse as sp

def get_options():
    parser = argparse.ArgumentParser(prog='parse_fragments.py')
    parser.add_argument('-i', '--input', help='BED file with fragments', required=True)
    parser.add_argument('-o', '--output', help='Output AnnData file', default='output.h5ad')
    parser.add_argument('-w', '--whitelist', help='List of barcodes')
    parser.add_argument('-b', '--bed_file', help='BED file with intervals for summarization')
    parser.add_argument('-g', '--genome_size', help='File with chromosome sizes')
    parser.add_argument('-D', '--count_duplicated', help='Count duplicated reads', action='store_true')
    parser.add_argument('-s', '--bin_size', help='Size of bins', type=int, default=5000)    
  
    options = parser.parse_args()
    return options

def regions_from_bed(bed_file):
    var_names = []
    regions = {}
    n_region = 0
    for line in open(bed_file):
        if line.startswith('#'):
            continue
        t = line.split()
        chrom, start, end = t[:3]
        if not chrom in regions:
            regions[chrom] = Intersecter()
        regions[chrom].add_interval(Interval(int(start), int(end), value=n_region))
        n_region += 1
        var_names.append(f"{chrom}:{start}-{end}")
    return (regions, var_names)        
        
def regions_from_bins(genome_size, bin_size=5000):
    var_names = []
    regions = {}
    n_region = 0
    for line in open(genome_size):
        if line.startswith('#'):
            continue
        t = line.split()
        chrom, chrom_size = t[:2]
        chrom_size = int(chrom_size)
        if not chrom in regions:
            regions[chrom] = Intersecter()
        for start in np.arange(0, chrom_size, bin_size):
            end = min(start + bin_size, chrom_size)
            regions[chrom].add_interval(Interval(int(start), int(end), value=n_region))
            n_region += 1
            var_names.append(f"{chrom}:{start}-{end}")
    return (regions, var_names)        

def parse_fragments():
    options = get_options()

    if not options.bed_file and not options.genome_size:
        sys.stderr.write("You should provide either a BED file or genome size\n")
        sys.exit(1)

    if options.bed_file:
        regions, var_names = regions_from_bed(options.bed_file)
    elif options.genome_size:
        regions, var_names = regions_from_bins(options.genome_size, options.bin_size)
    
    if options.whitelist:
        whitelist = [x.split()[0] for x in open(options.whitelist)]
    else:
        sys.stderr.write("You did not provided a whitelist, now extracting from BED file\n")
        sys.stderr.write("This is dangerous as it does not grant consistency\n")
        whietlist = set()
        for line in open(options.input):
            t = line.split()
            whitelist.add(t[3])
        whitelist = list(whitelist)
    bidx = dict([(whitelist[x], x) for x in range(len(whitelist))])
    
    # do the stuff
    M = sp.lil_matrix((len(whitelist), len(var_names)))
    for line in tqdm(open(options.input)):
        chrom, start, end, cell_bc, dup_count = line.split()
        if not cell_bc in whitelist:
            continue
        if not chrom in regions:
            continue
        count = 1
        if options.count_duplicated:
            count = int(dup_count)
        ovlp = regions[chrom].find(int(start), int(end))
        r_idx = bidx[cell_bc]
        for ovl in ovlp:
            c_idx = ovl.value
            M[r_idx, c_idx] = M[r_idx, c_idx] + count
    
    M = sp.csr_matrix(M)
    adata = ad.AnnData(M)
    adata.obs_names = whitelist
    adata.var_names = var_names
    
    adata.write(options.output)
        
    
    


if __name__ == '__main__':
    parse_fragments()

import pandas as pd
import numpy as np
import scipy.sparse as sp
import anndata as ad
from tqdm import tqdm
import yaml
import sys

config_file = sys.argv[1]
sample = sys.argv[2]

config = yaml.safe_load(open(config_file))
bins = pd.read_table(config['bed_file'], header=None)
chroms = np.unique(bins[0])
offsets = dict.fromkeys(chroms)
for c in chroms:
    offsets[c] = np.where(bins[0] == c)[0][0]
barcodes = config['barcodes']
whitelist = [x.split()[0] for x in open(f"{sample}_whitelist.tsv")]
bidx = dict([(whitelist[x], x) for x in range(len(whitelist))])
bin_size = 5000
M = dict([(x, sp.lil_matrix((len(whitelist), len(bins)))) for x in barcodes.keys()])

for tn in barcodes:
    for bc in barcodes[tn]:
        print(tn, bc)
        for line in tqdm(open(f"{sample}_BC_{bc}.bed")):
            chrom, start, end, name, dup_count = line.split()
            if not chrom in offsets:
                continue
            os = int(start) // bin_size
            oe = int(end) // bin_size
            dup_count = int(dup_count)
            if os != oe:
                c_idx = os + offsets[chrom]
                r_idx = bidx[name]
                M[tn][r_idx, c_idx] = M[tn][r_idx, c_idx] + 1
                c_idx = oe + offsets[chrom]
                r_idx = bidx[name]
                M[tn][r_idx, c_idx] = M[tn][r_idx, c_idx] + 1
            else:
                c_idx = os + offsets[chrom]
                r_idx = bidx[name]
                M[tn][r_idx, c_idx] = M[tn][r_idx, c_idx] + 1
        
for tn in M:
    M[tn] = sp.csr_matrix(M[tn])

adata = ad.AnnData(M['tn5'])
for tn in M.keys():
    adata.layers[tn] = M[tn]
adata.obs_names = whitelist
adata.var_names = [f"{x[0]}:{x[1]}-{x[2]}" for x in bins.values]

adata.write(f"{sample}.h5ad")

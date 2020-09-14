import scanpy as sc
import numpy as np
import pandas as pd
import sys

prefix=sys.argv[1]

adata = sc.read_mtx("%s.mtx" % prefix)
adata.obs.index = pd.read_csv('%s.barcodes.txt' % prefix, header=None)[0].values
adata.var.index = pd.read_csv('%s.genes.txt' % prefix, header=None)[0].values
adata.var.loc[:, 'n_cells'] = np.array(np.sum(adata.X, axis=0), dtype=int).ravel()
adata.obs.loc[:, 'n_regions'] = np.array(np.sum(adata.X, axis=1), dtype=int).ravel()

adata.write("%s.h5ad" % prefix)

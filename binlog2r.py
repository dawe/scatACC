import sys
import argparse 
import scipy.sparse as sp
import numpy as np
import pandas as pd
import scanpy as sc
import pyranges as pr
import pywt

def get_options():
	parser = argparse.ArgumentParser(prog='cna_analysis.py')
	parser.add_argument('-i', '--input-file', help='AnnData File with binned counts')
	parser.add_argument('-a', '--annotation', help='Bin annotation file')
	parser.add_argument('-p', '--prefix', help='Prefix for output file')
	parser.add_argument('--keep-bg', help='Keep background data (if any)', action='store_true')

	
	options = parser.parse_args()

	if options.bin_size > options.window_size:
		sys.stderr.write("Step size cannot be larger than window size\nSetting it to %d" % options.window_size)
		options.bin_size = options.window_size

	return options


def main():
	options = get_options()

	adata = sc.read(options.input_file)
	binner = pd.read_table(options.annotation, header=None)

	if not options.keep_bg:
		cells = [x for x in adata.obs_names if not x.startswith("Background")]
		adata = adata[cells]


	binner.columns = ['chromosome', 'start', 'end', 'gc_content', 'mapability', 'idx', 'idx_5k']
	
	binned_data = np.zeros((adata.shape[0], len(binner)))

    for x in range(len(binner)):
        V = adata[:, [int(x) for x in binner.iloc[x, -1].split(',')]].X.sum(axis=1).A1
        binned_data[:, x] = V

    binned_data = anndata.AnnData(binned_data)
    binned_data.obs_names = adata.obs_names
    for O in adata.obs.columns:
    	binned_data.obs[O] = adata.obs[O]
    	
    binned_data.var_names = [f"{x[0]}:{x[1]}-{x[2]}" for x in binner.values]
    nz = np.sum(binned_data.X != 0, axis=1)
    binned_data.obs['dropout'] = (binned_data.X.shape[1] - nz) / binned_data.X.shape[1] 
    binned_data.obs['n_reads'] = binned_data.X.sum(axis=1)

    binned_data.var['gc_content'] = binner['gc_content'].values
    binned_data.var['mapability'] = binner['mapability'].values

    binned_data.write(f"{options.prefix}_bin_counts.h5ad")
			
    binned_data.var.fillna(0, inplace=True)

    binned_data.X = np.log2( binned_data.X / np.mean(adata.X, axis=1)[:, None])

    binned_data.X[np.bitwise_not(np.isfinite(binned_data.X))] = -np.ptp(binned_data.X[np.isfinite(binned_data.X)])/2

    binned_data = binned_data.copy().T

    sc.pp.regress_out(binned_data, ['gc_content', 'mapability'])

    sc.pp.scale(test, max_value=4)

    binned_data = binned_data.copy().T
    binned_data.write(f"{options.prefix}_log2r.h5ad")


if __name__ == '__main__':
	main()


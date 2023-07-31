import scanpy as sc
import numpy as np
import pandas as pd
import sys
import anndata as ad
import argparse


def get_options():
    parser = argparse.ArgumentParser(prog='process_share.py')
    parser.add_argument('-s', '--spliced', help='mtx file with spliced counts')
    parser.add_argument('-u', '--unspliced', help='mtx file with unspliced counts')
    parser.add_argument('-v', '--var_names', help='file with var_names (genes)')
    parser.add_argument('-b', '--obs_names', help='file with obs_names (cells)')
    parser.add_argument('-p', '--prefix', help='Prefix for output files')
    parser.add_argument('-S', '--star_base', help='STAR out dir', default='')

    options = parser.parse_args()
  
    return options


def parse_mtx():

    options = get_options()
    if options.star_base:
        basedir = f'{options.star_base}/Gene/raw/'
        spliced_mtx = f'{basedir}/matrix.mtx'
        basedir = f'{options.star_base}/Velocyto/raw/'
        unspliced_mtx = f'{basedir}/matrix.mtx'
    else:
        spliced_mtx = options.spliced
        unspliced_mtx = options.unspliced       

    try:
        sdata = sc.read_mtx(spliced_mtx).T
    except FileNotFoundError:
        sys.stderr.write("Spliced counts file not found")
        sys.exit(1)

    try:
        udata = sc.read_mtx(unspliced_mtx).T
    except FileNotFoundError:
        sys.stderr.write("Unspliced counts file not found")
        if sdata:
           sys.stderr.write("Continuing anyway")
           udata = ''
        else:
            sys.exit(1)

    adata = ad.AnnData(sdata.X)
    if udata:
        adata.layers['spliced'] = sdata.X
        adata.layers['unspliced'] = udata.X
        
    if options.star_base:
        adata.obs_names = [x.strip() for x in open(f'{basedir}/barcodes.tsv')]
        import pandas as pd
        genes = pd.read_table(f'{basedir}/features.tsv', index_col=0, header=None)
        adata.var_names = list(genes.index)
        adata.var['Symbol'] = genes[1]
    else:    
        adata.obs_names = [x.strip() for x in open(options.obs_names)]
        adata.var_names = [x.strip() for x in open(options.var_names)]

    adata.var.loc[:, 'n_cells'] = np.array(np.sum(adata.X, axis=0), dtype=int).ravel()
    adata.obs.loc[:, 'n_regions'] = np.array(np.sum(adata.X, axis=1), dtype=int).ravel()

    adata.write(f'{options.prefix}.h5ad')


if __name__ == '__main__':
    parse_mtx()


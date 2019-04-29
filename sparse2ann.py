import numpy as np 
import scipy.sparse as sp
import anndata
import argparse
import sys

def get_options():
  parser = argparse.ArgumentParser(prog='sparse2ann.py')
  parser.add_argument('-p', '--peaks-file', help='Any BED-like file with peaks')
  parser.add_argument('-s', '--sparse', help='Sparse matrix file name (npz)')
  parser.add_argument('-a', '--anndata', help='Anndata file name (h5ad)')
  parser.add_argument('-i', '--inverse-transform', help='Perform inverse transformation (AnnData -> sparse)', action='store_true')
  parser.add_argument('-A', '--keep-all', help='Include all data, do not skip features with zero-coverage', action='store_true')
  
  options = parser.parse_args()
  if not options.sparse and not options.anndata:
    # for the time being, default output is bed
    sys.stderr.write("You should provide at least one file as input\n")
    sys.exit(1)

  if not options.anndata and (options.sparse and not options.inverse_transform):
  	sys.stderr.write("If you intend to perform transformation without any peak list, output will lack any annotation\n")

  if not options.sparse:
    option.sparse = 'output.npz'

  if not options.anndata:
    options.anndata = 'output.h5ad' 	

  return options


def main():
  options = get_options()

  if not options.inverse_transform:

    loader = np.load(options.sparse)
    try:
      count_matrix = sp.csr_matrix((loader['data'], loader['indices'], loader['indptr']), shape=loader['shape'])
    except ValueError:
      count_matrix = sp.csc_matrix((loader['data'], loader['indices'], loader['indptr']), shape=loader['shape'])
    bc_list = loader['bc_list']

    if count_matrix.shape[0] > count_matrix[1]:
      # we probably have data in regions by cells, we need it the other way
      count_matrix = sp.csr_matrix(count_matrix.T) #convert

    if options.keep_all:
      mask = np.ones(count_matrix.shape[1], dtype=bool)
    else:
      mask = np.array(np.sum(count_matrix, axis=0) > 0).ravel()

    
    if options.peaks_file:
      regions = []
      for line in options.peaks_file:
        t = line.split()	
        chrom, start, end = t.split()[:3]
        r_id = "%s:%s-%s" % (chrom, start, end)
        regions.append(r_id)
      regions = np.array(regions)  
    else:
      regions = np.arange(count_matrix.shape[1])  
    
    regions = regions[mask]
    n_cells = pd.DataFrame(np.array(np.sum(count_matrix > 0, axis=0)).ravel(), index=regions, columns=['n_cells'])
    n_regions = pd.DataFrame(np.array(np.sum(count_matrix > 0, axis=1)).ravel(), index=bc_list, columns=['n_regions'])
    adata = anndata.AnnData(count_matrix, obs=n_regions, var=n_cells)
    adata.write(options.anndata)

  else:
  	adata = anndata.read(options.anndata)
    count_matrix = 	adata.X
    bc_list = np.array(adata.obs.index)
    np.savez(options.sparse, data = count_matrix.data, indices=count_matrix.indices, 
             indptr=count_matrix.indptr, shape=count_matrix.shape, bc_list=bc_list)
    



if __name__ == '__main__':
  main()


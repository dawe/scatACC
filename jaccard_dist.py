import sys
import numba
import numpy as np
import scipy.sparse as sp
import anndata

@numba.njit(parallel=True, debug=False)
def jaccard_dist(indptr,indices):
  n = len(indptr - 1)
  out = np.zeros((n, n))
  for i in numba.prange(n - 1):
    A = set([indices[x] for x in numba.prange(indptr[i], indptr[i + 1])])
    for j in numba.prange(i + 1, n):
      B = set([indices[x] for x in numba.prange(indptr[j], indptr[j + 1])])
      I = 1. * len(A.intersection(B))
      U = 1. * len(A.union(B))
      if U == 0:
        J = 0.
      else:
        J = I / U
      out[i, j] = 1 - J   
      out[j, i] = out[i, j]  
  return out      



adata = anndata.read_h5ad(sys.argv[1])
count_matrix = adata.X

dm = jaccard_dist(count_matrix.indptr, count_matrix.indices)

np.savez('dm.npz', dm)
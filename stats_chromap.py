import sys
import os
import gzip
import scipy.optimize
import numpy as np
from scipy.special import lambertw as W

def fs(c, n):
    return np.real(c*n / (c*W(-np.exp(-n/c)*n/c) + n))

bedfile = sys.argv[1]
outfile = sys.argv[2]

if bedfile.endswith('.bed.gz'):
    prefix = os.path.basename(bedfile).replace('.bed.gz', '')
else:
    prefix = os.path.basename(bedfile).replace('.bed', '')

print(f"Sample:\t{prefix}")

with open(outfile) as f:
    lines = f.read().split('\n')
for line in lines[-15:-2]:
    if line.startswith('Number'):
        if line[-1] == '.':
            line = line[:-1]
        print(line.replace(': ', ':\t'))

if bedfile.endswith('.bed.gz'):
    with gzip.open(bedfile) as f:
        dups = [int(x.decode('ascii').split()[-1]) for x in f]
else:
    with open(bedfile) as f:
        dups = [int(x.split()[-1]) for x in f]
u, n = np.unique(dups, return_counts=True)
N = np.sum(u*n)
C = np.sum(n)

E = np.max([0, int(fs(C, N))])
dr = 1 - C / N
sat = np.min([1, C / E])
print(f"Duplication Rate:\t{dr:.3f}")
print(f"Estimated Complexity:\t{E}")
print(f"Saturation:\t{sat:.3f}")

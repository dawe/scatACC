import warnings
warnings.filterwarnings('ignore')
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import scipy.signal
import os
import sys

plt.rcParams['font.size'] = 8

bedfile=sys.argv[1]
sample=bedfile.replace('.bed', '')

isizes = []
for line in open(bedfile):
    t = line.split()
    if t[0]=='chrM':
        continue
    isizes.append(int(t[2]) - int(t[1]))
isizes = np.array(isizes)

fig = plt.figure(figsize=(8, 3))
ax=plt.subplot(1, 2, 1)
nf = 10
bins=np.arange(1, 2000, nf)
H = plt.hist(isizes, bins=bins, density=True)
a = np.fft.rfft(H[0])
NS = len(H[0])
NFRE = np.sum(H[0][:12]) / np.sum(H[0][12:24])
xr = np.fft.rfftfreq(NS, 1/nf)
p = np.linalg.norm(np.vstack([np.real(a), np.imag(a)]), axis=0)
v = scipy.signal.find_peaks(p, prominence=0.1)
o = np.zeros_like(a)
w = np.arange(15)
o[w] = a[w]

S = np.fft.irfft(o)
S = S-S.min()
ax.plot(bins[:len(S)], S, c='C1')
ax.set_yticks([])
ax.set_ylim(0, 2*S.max())
ax.set_xlim(1, 1001)
Xt = [100, 300, 500, 800, 1000]
ax.set_xticks(Xt)
ax.set_xticklabels(Xt)
ax.set_title('Size Distribution')
ax = plt.subplot(1, 2, 2)
ax.plot(xr[:40], p[:40])
Xt = np.array([50, 75, 150, 300, 2000])
ax.set_xticks(100/Xt)
ax.set_xticklabels(Xt, rotation=90)
ax.set_title('Spectrum')
fig.suptitle(f'{sample} (NFR Enrichment: {NFRE:.3f})')
#fig.tight_layout()    
plt.savefig(f'{sample}.png', dpi=120, bbox_inches='tight')

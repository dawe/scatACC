import warnings
#warnings.filterwarnings('ignore')
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
import seaborn as sns
import sklearn.preprocessing
import natsort
import sys
plt.rcParams['image.cmap'] = 'coolwarm'

tn_barcodes = {'tn5':"""CGTACTAG
TCCTGAGC
TCATGAGC
CCTGAGAT""".split(), 'tnH':"""TAAGGCGA
GCTACGCT
AGGCTCCG
CTGCGCAT""".split()}

df = pd.read_table(sys.argv[1], header=None)

#samples = df.columns[4:]
samples = [x.strip() for x in open(sys.argv[2])]
df.columns=['chrom', 'start', 'end', 'state'] + samples
p = df.groupby('state')[samples].agg('sum')
samples_sum = list(set([x.split('/')[1] for x in samples]))
samples_sum.sort()
samples = [f'{s}_{t}' for s in samples_sum for t in ['tn5', 'tnH']]

df2 = pd.DataFrame(0, index=df.index, columns=samples_sum)

for s in samples_sum:
    sdf = df.filter(like=f'{s}_BC_')
    for t in ['tn5', 'tnH']:
        tn_mask = np.where([len(set(x.split('_')).intersection(tn_barcodes[t])) > 0 for x in sdf.columns])[0]
        df2[f'{s}_{t}'] = sdf.iloc[:, tn_mask].sum(1)

df2['state'] = df['state']
p = df2.groupby('state')[samples].agg('sum')


#states = ['1_TssA', '2_TssAFlnk', '3_TxFlnk', '4_Tx', '5_TxWk', '6_EnhG', '7_Enh', '8_ZNF/Rpts', '9_Het', '10_TssBiv', '11_BivFlnk', '12_EnhBiv', '13_ReprPC', '14_ReprPCWk', '15_Quies']

#p.index=states
p = p / p.sum()

df = p.iloc[:, 1::2] / p.iloc[:, ::2].values
desc = df.T.describe().T
desc.to_csv("State_Enrichment_Stats.txt", sep="\t")
df.to_csv("State_Enrichment.txt", sep="\t")
m = desc['mean']
s = desc['std']

idx = m.sort_values().index
m = m[idx]
s = s[idx]
fig = plt.figure()
ax = fig.add_subplot(1,1,1)

sns.boxplot(data = df.T[idx], orient='h', palette=plt.cm.hsv([int(x.split('_')[0])/15 for x in df.T[idx].columns]), fliersize=0, ax=ax)
sns.swarmplot(data = df.T[idx], orient='h', color='black', ax=ax)


#yticks(range(len(df)), m.index)
ax.set_title("State Enrichment")
ax.vlines(1, -1, 15)
ax.set_ylim(-.5, 14.5)
ax.set_xlabel("tnH / tn5")
plt.tight_layout()
fig.savefig("State_Enrichment.png", dpi=300)

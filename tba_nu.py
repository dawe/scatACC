import sys
import numpy as np
import numba
#from joblib import Parallel, delayed
#import multiprocessing

# tba.py fasta_file tf.meme 

num_cores = 4


def reverse_complement(s):
    ra = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
    return ''.join([ra[x] for x in s[::-1]])

tf_d = {}
for line in open(sys.argv[2]):
    if line.startswith('MOTIF'):
        t = line.split()
        tf_name = t[1]
    if line.startswith('letter-probability'):
        t = line.split()
        tf_len = int(t[5])
        tf_d[tf_name] = np.zeros((tf_len, 4))
        lc = 0
    if line.startswith('0') or line.startswith('1'):
        tf_d[tf_name][lc] = [float(x) for x in line.split()]
        lc += 1

for tf_name in tf_d:
    tf_d[tf_name] = np.log(tf_d[tf_name] + 1e-9)    

bg_p = {'A':2.523e-01,
'C': 2.477e-01,
'G': 2.477e-01,
'T': 2.523e-01,
}
nt_idx = {'A':0, 'C':1, 'G':2, 'T':3}
for nt in 'ACGT':
    bg_p[nt] = np.log(bg_p[nt])
bg_pv = np.array(list(bg_p.values()))    

#@numba.jit(numba.float64(numba.int8[:], numba.float64[:, :], numba.float64[:]))
#@numba.jit(numba.float64(numba.int32[:], numba.float64[:, :], numba.float64[:]))
@numba.jit(nopython=True)
def n_affinity(seq, tf_mat, bg_pv):
    p = np.zeros(len(seq))
    for x in range(len(tf_mat)):
        p[x] = tf_mat[x, seq[x]]
#    l_tf = len(tf_mat)
#    p = tf_mat[np.arange(l_tf), seq]
    bg = bg_pv[seq]
    return np.sum(p - bg)

#@numba.jit(numba.float64(numba.int32[:, :],numba.int32[:, :],numba.float64[:, :],numba.float64[:], numba.int8), nopython=True)
@numba.jit(nopython=True)
def n_tba(seq_array, rseq_array, tf_mat, bg_pv, log=True):
#    aff_fwd = np.zeros(len(seq_array))
#    aff_rev = np.zeros(len(seq_array))
    tba_arr = np.zeros(len(seq_array))
    for x in range(len(tba_arr)):
        aff_fwd = n_affinity(seq_array[x], tf_mat, bg_pv)
        aff_rev = n_affinity(rseq_array[x], tf_mat, bg_pv)
        tba_arr[x] = np.exp(max(aff_fwd, aff_rev))
    tba = np.sum(tba_arr)
    if log:
        return np.log(tba)
    return tba    

def replace_random(seq):
    alphabet = 'ACGT'
    nts = [x for x in seq]
    for x in range(len(seq)):
        if nts[x] == 'N':
            nts[x] = alphabet[np.random.randint(0, 4)]
    return ''.join(nts)    

tba_seq = {}
tf_names = list(tf_d.keys())
import time

all_seqs = {}
for line in open(sys.argv[1]):
    if line.startswith('>'):
        s_id = line.strip()[1:]
    else:
        seq = line.strip().upper()
        if 'N' in seq:
            seq = replace_random(seq)
        all_seqs[s_id] = seq

print(f'read {len(all_seqs)} sequences')

tba_seq = {}
nq = 0
for ns, s_id in enumerate(all_seqs):
    seq_name = s_id.split('::')[0] #built like this >merged_500bp::summit
    if not seq_name in tba_seq:
        S = time.time()
        tba_seq[seq_name] = np.zeros(len(tf_names))
        nq += 1
    seq = all_seqs[s_id]
    seq_idx = np.array([nt_idx[x] for x in seq])
    rseq_idx = np.array([nt_idx[x] for x in reverse_complement(seq)])
    tba = np.zeros(len(tf_names))
    l_seq = len(seq)
    for tn, tx in enumerate(tf_names):
        l_tf = len(tf_d[tx])
        seq_array = np.array([seq_idx[s:s + l_tf] for s in range(l_seq - l_tf)])
        rseq_array = np.array([rseq_idx[s:s + l_tf] for s in range(l_seq - l_tf)])    
        tba[tn] = n_tba(seq_array, rseq_array, tf_d[tx], bg_pv, log=False) 
    tba_seq[seq_name] += tba
    if ns % 100 == 0:
        print(nq, ns, time.time() - S)

for seq_name in tba_seq:
    tba_seq[seq_name] = np.log(tba_seq[seq_name])
    
fout = sys.argv[1].replace('.fa', '.pickle')  
import pandas as pd          
pd.DataFrame(tba_seq, index=tf_names).to_pickle(fout)

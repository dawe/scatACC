import sys
import numpy as np
import numba
from tqdm import tqdm
import argparse
import pysam
import pandas as pd          

#from joblib import Parallel, delayed
#import multiprocessing

# tba.py fasta_file tf.meme 

num_cores = 12

def get_options():
    parser = argparse.ArgumentParser(prog='mmc.py')
    parser.add_argument('-f', '--fasta', help='Fasta file with sequences', required=True)
    parser.add_argument('-m', '--matrix', help='MEME formatted file with matrices', required=True)
    parser.add_argument('-o', '--output', help='Output file name')
    parser.add_argument('-p', '--pad', help='Pad sequences to this length', default=0)
    parser.add_argument('-b', '--background', help='Tab separated file with nucleoutide probability')
    
  
    options = parser.parse_args()
  
    return options


def reverse_complement(s):
    ra = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
    return ''.join([ra[x] for x in s[::-1]])


@numba.jit(numba.float64(numba.int8[:], numba.float64[:, :], numba.float64[:]), nopython=True)
#@numba.jit(nopython=True)
def n_affinity(seq, tf_mat, bg_pv):
    p = np.zeros(seq.size)
    for x in range(seq.size):
        p[x] = tf_mat[x, seq[x]]
    bg = bg_pv[seq]
    return np.sum(p - bg)

@numba.jit(numba.float64(numba.int8[:,:],numba.int8[:,:],numba.float64[:,:],numba.float64[:], numba.boolean), nopython=True)
#@numba.jit(nopython=True)
def n_tba(seq_array, rseq_array, tf_mat, bg_pv, log=True):
    tba_arr = np.zeros(seq_array.shape[0])
    for x in range(tba_arr.size):
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

def read_matrix(file):
    tf_d = {}
    for line in open(file):
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

    return tf_d


def main():
    options = get_options()

    tf_d = read_matrix(options.matrix)

    bg_p = {}
    if not options.background:
        bg_p = {'A':2.5e-01,
                'C':2.5e-01,
                'G':2.5e-01,
                'T':2.5e-01,
               }
    else:
        for line in open(options.background):
            t = line.split()
            bg_p[t[0]] = float(t[1])

    nt_idx = {'A':0, 'C':1, 'G':2, 'T':3}
    for nt in 'ACGT':
        bg_p[nt] = np.log(bg_p[nt])
    bg_pv = np.array(list(bg_p.values()))    


    tba_seq = {}
    tf_names = list(tf_d.keys())
    n_tf = len(tf_names)

    seq_fh = pysam.FastxFile(options.fasta)
    for entry in tqdm(seq_fh, unit='sequence'):
        seq = entry.sequence.upper()
        seq_name = entry.name
        if not seq_name in tba_seq:
            # initialize
            tba_seq[seq_name] = np.zeros(n_tf)
        if options.pad > 0:
            l = len(seq)
            p = options.pad - l
            seq = seq + ('N' * p)
        if 'N' in seq:
            seq = replace_random(seq)
        
        seq_idx = np.array([nt_idx[x] for x in seq])
        rseq_idx = np.array([nt_idx[x] for x in reverse_complement(seq)])
        tba = np.zeros(len(tf_names))
        l_seq = len(seq)
        for tn, tx in enumerate(tf_names):
            l_tf = len(tf_d[tx])
            seq_array = np.array([seq_idx[s:s + l_tf] for s in range(l_seq - l_tf)], dtype=np.int8)
            rseq_array = np.array([rseq_idx[s:s + l_tf] for s in range(l_seq - l_tf)], dtype=np.int8)    
            tba[tn] = n_tba(seq_array, rseq_array, tf_d[tx], bg_pv, log=False) 
        tba_seq[seq_name] += tba

    for seq_name in tba_seq:
        # log transform
        tba_seq[seq_name] = np.log(tba_seq[seq_name])

    if not options.output:
        fout = options.fasta.replace('.fa', '.pickle')        
    pd.DataFrame(tba_seq, index=tf_names).T.to_pickle(fout)

if __name__ == '__main__':
  main()

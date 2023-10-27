import HTSeq
import sys
import argparse
from tqdm import tqdm

def get_options():
    parser = argparse.ArgumentParser(prog='hash_match.py')
    parser.add_argument('-o', '--output', help='Name of output file', default='hash_match.csv')
    parser.add_argument('-1', '--read1', help='Fastq file with R1', required=True)
    parser.add_argument('-2', '--read2', help='Fastq file with R2', required=True)
    parser.add_argument('-W', '--whitelist_rna', help='UMI-tools whitelist for RNA', required=True)
    parser.add_argument('-H', '--whitelist_hash', help='UMI-tools whitelist for Hash', required=True)
    parser.add_argument('-n', '--n_seq', help='Max number of sequences to process (for debugging)', default=0, type=int)
    
    options = parser.parse_args()
    
    return options


def process_tables():
    options = get_options()
    
    hash_dict = {}
    bc_dict = {}
    
    hmatch = {}
    
    for line in open(options.whitelist_rna):
        base, alt, counts, _ = line.split('\t')
        base = bytes(base, encoding='ascii')
        bc_dict[base] = base
        for a in alt.split(','):
            a = bytes(a, encoding='ascii')
            bc_dict[a] = base
    
    for line in open(options.whitelist_hash):
        base, alt, counts, _ = line.split('\t')
        base = bytes(base, encoding='ascii')
        hash_dict[base] = base
        for a in alt.split(','):
            a = bytes(a, encoding='ascii')
            hash_dict[a] = base
    
    
    
    r1 = HTSeq.FastqReader(options.read1)
    r2 = HTSeq.FastqReader(options.read2)    
    
    read_iterator = zip(r1, r2)

    n_tot = 0    
    for item in tqdm(read_iterator):
        if options.n_seq > 0 and n_tot == options.n_seq:
            break
    
        rna_bc = item[0].seq[:16]
        cell_hash = item[1].seq[:8]
        
        if rna_bc in bc_dict:
            rna_bc = bc_dict[rna_bc]
        if cell_hash in hash_dict:
            cell_hash = hash_dict[cell_hash]
        
        if not rna_bc in hmatch:
            hmatch[rna_bc] = {cell_hash:1}
        elif not cell_hash in hmatch[rna_bc]:
            hmatch[rna_bc][cell_hash] = 1
        else:
            hmatch[rna_bc][cell_hash] += 1
        
        n_tot += 1
         
    with open(options.output, 'w') as fout:
        for bc in hmatch:
            for ha in hmatch[bc]:
                fout.write(f'{bc.decode()}\t{ha.decode()}\t{hmatch[bc][ha]}\n')        
        




if __name__ == '__main__':
    process_tables()

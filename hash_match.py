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
    parser.add_argument('-k', '--skip_unmatched', help='Skip entries not listed in any whitelist', action='store_true')
    parser.add_argument('-D', '--skip_duplicated', help='Skip entries with duplicated UMI', action='store_true')
    parser.add_argument('-u', '--umi_length', help='Length of the UMI after the CB', default=12, type=int)
    
    options = parser.parse_args()
    
    return options

def hamming(a, b):
    la = len(a)
    lb = len(b)
    if la != lb:
        return la
    return sum([a[x] != b[x] for x in range(len(a))])


def process_tables():
    options = get_options()
    
    hash_dict = {}
    bc_dict = {}
    
    hmatch = {}
    
    for line in open(options.whitelist_rna):
        t = line.split('\t')
        base = bytes(t[0], encoding='ascii')
        bc_dict[base] = base
        for a in t[1].split(','):
            a = bytes(a, encoding='ascii')
            bc_dict[a] = base
    
    for line in open(options.whitelist_hash):
        t = line.split('\t')
        base = bytes(t[0], encoding='ascii')
        hash_dict[base] = base
        for a in t[1].split(','):
            a = bytes(a, encoding='ascii')
            hash_dict[a] = base

   # apparently umitools has issues with short barcodes and does not collapse
   # correctly these 8 bp        
    
    
    
    r1 = HTSeq.FastqReader(options.read1)
    r2 = HTSeq.FastqReader(options.read2)    
    
    read_iterator = zip(r1, r2)

    n_tot = 0    
    umi_start = 16
    umi_end = umi_start + options.umi_length
    umi_counter = {}
    for bc in bc_dict.keys():
        umi_counter[bc] = set()
    for item in tqdm(read_iterator):
        if options.n_seq > 0 and n_tot == options.n_seq:
            break
        n_tot += 1
         
    
        rna_bc = item[0].seq[:16]
        rna_umi = item[0].seq[umi_start:umi_end]
        cell_hash = item[1].seq[:8]
        
        skip = 0 
        
        if rna_bc in bc_dict:
            rna_bc = bc_dict[rna_bc]
        else:
            skip += 1
        if cell_hash in hash_dict:
            cell_hash = hash_dict[cell_hash]
        else:
            skip += 1
    
        if skip > 0 and options.skip_unmatched:
            continue
        
        if rna_umi in umi_counter[rna_bc] and options.skip_duplicated:
            continue
        else:
            umi_counter[rna_bc].add(rna_umi)


        if not rna_bc in hmatch:
            hmatch[rna_bc] = {cell_hash:1}
        elif not cell_hash in hmatch[rna_bc]:
            hmatch[rna_bc][cell_hash] = 1
        else:
            hmatch[rna_bc][cell_hash] += 1
   


        
    with open(options.output, 'w') as fout:
        fout.write("cell_barcode\thash\tcount\n")
        for bc in hmatch:
            for ha in hmatch[bc]:
                fout.write(f'{bc.decode()}\t{ha.decode()}\t{hmatch[bc][ha]}\n')        
        




if __name__ == '__main__':
    process_tables()

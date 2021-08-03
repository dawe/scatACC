import sys
import argparse
import editdistance as ed
import pandas as pd
import tqdm

def get_options():
  parser = argparse.ArgumentParser(prog='mmc.py')
  parser.add_argument('-w', '--whitelist', help='Whitelist file from UMI_tools', required=True)
  parser.add_argument('-b', '--barcodes', help='Precompiled barcode list', required=True)
  parser.add_argument('-o', '--output', help='Output whitelist file')
  parser.add_argument('-r', '--reverse_complement', help='Whitelist is in reversec complement', action='store_true')
  parser.add_argument('-d', '--distance', help='Max edit distance allowed', default=1)  
  parser.add_argument('-E', '--exhaustive_search', help='Search for matches that are not perfect', action='store_true')
  
  options = parser.parse_args()
  
  return options


def reverse_complement(s):
  ab = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', ',':','}
  return ''.join([ab[x] for x in s[::-1]])

def rc_df(wl_df):
  # reverse complement all entries
  rc_df = pd.DataFrame(wl_df.values, index=[reverse_complement(x) for x in wl_df.index], columns=wl_df.columns)
  # reverse complement everything, doesn't matter the order
  rc_df[1] = [reverse_complement(x) for x in wl_df[1]]
  # just reverse the counts
  rc_df[3] = [','.join(x.split(',')[::-1]) for x in wl_df[3]]
  return rc_df
  

def reshuffle(new_bc, old_bc, synonyms, t_counts, syn_counts):
  count_list = syn_counts.split(',')
  syn_sum = sum([int(x) for x in count_list])
  bc_count = t_counts - syn_sum
  syn_l = synonyms.split(',')
  syn_idx = syn_l.index(new_bc)
  syn_l[syn_idx] = old_bc
  count_list[syn_idx] = str(bc_count)
  synonyms = ','.join(syn_l)
  syn_counts = ','.join(count_list)
  return [new_bc, synonyms, t_counts, syn_counts]

def reassign(new_bc, old_bc, synonyms, t_counts, syn_counts):
  syn_sum = sum([int(x) for x in syn_counts.split(',')])
  bc_count = t_counts - syn_sum
  synonyms = f'{old_bc},{synonyms}'
  syn_counts = f'{bc_count},{syn_counts}'
  return [new_bc, synonyms, t_counts, syn_counts]
  


def main():
  options = get_options()
  if options.output:
    fout = options.output 
  else:
    fout = sys.stdout
    
  wl_df = pd.read_table(options.whitelist, header=None, index_col=0) # from UMI_tools only...
  bc_list = set([x.strip() for x in open(options.barcodes)]) # from 10x?
  
  # convert to rc if needed
  if options.reverse_complement:
    wl_df = rc_df(wl_df)
  # first of all we keep exact matches in the index, to reduce the checks
  
  known = bc_list.intersection(wl_df.index)
  unknown = set(wl_df.index) - bc_list
  
  sys.stderr.write(f'{len(known)}/{len(wl_df)} perfect barcode match\n')
  
  sys.stderr.write(f'Analyzying the remaining {len(unknown)} barcodes\n')
  reassign_idx = []
  reassign_val = []
  for bc in tqdm.tqdm(unknown):
    synonyms, t_counts, syn_counts = wl_df.loc[bc]
    syn_s = set(synonyms.split(','))
    candidates = list(syn_s.intersection(bc_list))
    if len(candidates) > 0:
      # reshuffle
      _p = reshuffle(candidates[0], bc, synonyms, t_counts, syn_counts)
      reassign_idx.append(_p[0])
      reassign_val.append(_p[1:])
    elif options.exhaustive_search:
      # try harder
      for bc_s in syn_s:
        candidates = [x for x in bc_list if ed.eval(bc_s, x) <= options.distance] #can be any of these, actually
        if len(candidates) > 0:
          _p = reassign(candidates[0], bc, synonyms, t_counts, syn_counts)
          reassign_idx.append(_p[0])
          reassign_val.append(_p[1:])
          break
  
  if len(reassign_idx) > 0:
    reassigned = pd.DataFrame(reassign_val, index=reassign_idx, columns=wl_df.columns)
    out_df = pd.concat([reassigned, wl_df.loc[known]], axis=0)
  else:
    out_df = wl_df.loc[known]
  if options.reverse_complement:
    # if reverse_complement, then re-reverse otherwise BC cannot be found in
    # reads
    out_df = rc_df(out_df)
  out_df.to_csv(fout, sep="\t", header=None)
  
  
  


if __name__ == '__main__':
  main()
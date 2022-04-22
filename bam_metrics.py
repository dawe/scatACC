import os
import pysam
import sys

fh = pysam.AlignmentFile(sys.argv[1], 'rb')

mapped = 0
dups = 0
total = 0
bg_read = 0

for r in fh:
	if r.is_secondary or r.is_supplementary:
		continue
	total += 1
	if r.is_duplicate:
		dups += 1
	if not r.is_unmapped:
		mapped += 1
	if dict(r.tags)['RG'].startswith("Back"):
		bg_read += 1

fn=os.path.basename(sys.argv[1])
print(f'File\t{fn}')
print(f'Total\t{total}')
print(f'Mapped\t{mapped}')
print(f'Duplicated\t{dups}')
print(f'In Cell\t{total - bg_read}')
r=mapped/total
print(f'%Mapped\t{r:.4f}')
r=dups/mapped
print(f'%Duplicated\t{r:.4f}')
r=1 - (bg_read / total)
print(f'%In Cell\t{r:.4f}')

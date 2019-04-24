import sys
import pysam


wl_file = sys.argv[1]
peak_file = sys.argv[2]
bam_file = sys.argv[3]


#read whitelist
counter = {}
for line in open(wl_file):
  t = line.split()
  counter[t[0]] = 0
bc_list = list(counter.keys())

spool = ['chrom', 'start', 'end'] + bc_list
sys.stdout.write("\t".join(spool) + "\n")

fin = pysam.Samfile(bam_file, 'rb')
for line in open(peak_file):
  if line.startswith('#'):
    continue
  chr, start, end = line.split()[:3]
  for alignment in fin.fetch(chr, int(start), int(end)):
    if alignment.is_proper_pair and alignment.mapq >= 15 and alignment.is_read1:
      counter[dict(alignment.tags)['CB']] += 1
  spool = [chr, start, end] + ["%d" % counter[x] for x in bc_list]
  sys.stdout.write("\t".join(spool) + "\n")
  counter = dict([(x,0) for x in bc_list])

fin.close()    
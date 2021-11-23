import sys
import pysam

if len(sys.argv) < 3:
    print("You should specify BAM in and BAM out")
    sys.exit(1)

fin = sys.argv[1]
fout = sys.argv[2]

bam_in = pysam.AlignmentFile(fin, 'rb')
in_header = bam_in.header

bam_out = pysam.AlignmentFile(fout, 'wb', header=in_header)
rg_sm = {}

for record in in_header['RG']:
    rg_sm[record['ID']] = record['SM']

for record in bam_in:
    tags = record.tags.copy()
    rg = [x[1] for x in tags if x[0]=='RG'][0]
    sm = rg_sm[rg]
    if sm != 'Background':
	    tags.append(('CB', rg_sm[rg]))
	    record.tags = tags
    bam_out.write(record)

bam_in.close()
bam_out.close()    

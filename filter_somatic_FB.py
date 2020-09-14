import sys
import argparse
import cyvcf2
import pandas as pd


def _check_lods(record, tum_idx, norm_idx, tumor_thresh=35, normal_thresh=35):
    """Ensure likelihoods for tumor and normal pass thresholds.
    Skipped if no FreeBayes GL annotations available.
    """
    try:
        tumor_gls = record.gt_phred_ll_homref[tum_idx], record.gt_phred_ll_het[tum_idx], record.gt_phred_ll_homalt[tum_idx]
        tumor_lod = min(tumor_gls[i] - tumor_gls[0] for i in range(1, len(tumor_gls)))
      # No GL information, no tumor call (so fail it)
    except IndexError:
        tumor_lod = -1.0
    try:
        normal_gls = record.gt_phred_ll_homref[norm_idx], record.gt_phred_ll_het[norm_idx], record.gt_phred_ll_homalt[norm_idx]
        normal_lod = min(normal_gls[i] - normal_gls[0] for i in range(1, len(normal_gls)))
    # No GL inofmration, no normal call (so pass it)
    except IndexError:
        normal_lod = normal_thresh
    return normal_lod >= normal_thresh and tumor_lod >= tumor_thresh

def _check_freqs(record, tum_idx, norm_idx, thresh_ratio=2.7):
    """Ensure frequency of tumor to normal passes a reasonable threshold.
    Avoids calling low frequency tumors also present at low frequency in normals,
    which indicates a contamination or persistent error.
    """
    tumor_freq = record.gt_alt_freqs[tum_idx]
    normal_freq = record.gt_alt_freqs[norm_idx]
    if tumor_freq == -1 or normal_freq == -1:
        return False
    return normal_freq <= 0.001 or normal_freq <= tumor_freq / thresh_ratio

def call_somatic(record, tum_idx, norm_idx, tumor_thresh = 35, normal_thresh = 35):
    # Thresholds are like phred scores, so 3.5 = phred35
    if record.gt_depths[tum_idx] == -1 or record.gt_depths[norm_idx] == -1:
        return False
    if _check_lods(record, tum_idx, norm_idx, tumor_thresh, normal_thresh) and _check_freqs(record, tum_idx, norm_idx):
        return True
    return False    
    
def main():
  option_parser = argparse.ArgumentParser(
  description="Based on Brad Chapman code",
  prog="filter_somatic_FB.py",
  epilog="For any question, write to cittaro.davide@hsr.it")
  option_parser.add_argument("--version", action="version", version="%(prog)s 0.1")
  option_parser.add_argument("--vcf", help="Multisample VCF file", action='store', type=str, required=True)
#  option_parser.add_argument("--name", help="name", action='store', default='script')

  options = option_parser.parse_args()  
  #opzioni
  
  if not options.vcf:
    vcf_parser = cyvcf2.VCF(sys.stdin)
  else:
    vcf_parser = cyvcf2.VCF(options.vcf)  

#  vcf_out = vcf.VCFWriter(sys.stdout, template=vcf_parser)  

  samples = vcf_parser.samples
  df = []
  index = []
  for record in vcf_parser:
    somatic = [False for x in samples]
    for X in range(len(samples)):
      for Y in range(len(samples)):
        if X != Y and call_somatic(record, X, Y):
          somatic[X] = True
    if sum(somatic) > 0:            
      df.append(somatic)
      index.append("%s:%d:%s:%s" % (record.CHROM, record.POS, record.REF, ','.join(record.ALT)))


  
  for record in vcf_parser:
    if  call_somatic(record, tum_idx, norm_idx):
      vcf_out.write_record(record)
    
if __name__ == '__main__':
  main()

import anndata as ad
import sys
import argparse


def get_options():
  parser = argparse.ArgumentParser(prog='sum_barcodes.py')
  parser.add_argument('-s', '--sample_name', help='Sample name', required=True)
  parser.add_argument('-i', '--input_files', help='Count files to be considered', nargs='+')
  parser.add_argument('-b', '--barcodes', help='Non-default barcode table')
  parser.add_argument('-X', '--Xdata', help='Layer to be used as anndata.X', default='tn5')  
  
  options = parser.parse_args()
  
  return options


def main():
  options = get_options()
  
  barcodes = {'CGTACTAG':'tn5',
              'TCCTGAGC':'tn5',
              'TCATGAGC':'tn5',
              'CCTGAGAT':'tn5',
              'TAAGGCGA':'tnH',
              'GCTACGCT':'tnH',
              'AGGCTCCG':'tnH',
              'CTGCGCAT':'tnH'
             }

  if options.barcodes:
    barcodes = {}
    for line in open(options.barcodes):
      t = line.split()
      barcodes[t[0]] = t[1]
      
  ad_l = dict.fromkeys(barcodes.values())
  for l in ad_l:
    bc_l = [x for x in barcodes if barcodes[x] == l]
    layer_files = [y for x in bc_l for y in options.input_files if x in y]
    ad_tmp = ad.read(layer_files[0])
    for f in layer_files[1:]:
      _X = ad.read(f)
      _X = _X[ad_tmp.obs_names]
      ad_tmp.X = ad_tmp.X + _X.X
    ad_l[l] = ad_tmp.copy()
  adata = ad_l[options.Xdata].copy()
  for l in ad_l:
    adata.layers[l] = ad_l[l].X

  adata.write(f'{options.sample_name}.h5ad')      


if __name__ == '__main__':
  main()
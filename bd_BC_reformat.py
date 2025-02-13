import pandas as pd
import sys
import numpy as np
import bgzip
#import editdistance as ed
from tqdm import tqdm
import argparse
import HTSeq


# struct R1
# NNNNNNNNNNN…
# struct R2
# CAAGCGTTGGCTTCTCGCATCT AGTGGTCA ATCCACGTGCTTGAG AGGCCAGAGCATTCG AACGCTTA
# ---------------------- 22222222 --------------- --------------- 11111111

# in RNA
# GAAGCGTTGGCTTCTCGCATCT CAACCACA ATCCACGTGCTTGAG AGGCCAGAGCATTCG ACATTGGC GTGGCCGATGTTTCGCATCGGCGTACGA CTTAGTGGGT ATTTTTTTTTTTTTTTGTTTATGGGGTTTTTTTTGGTTTTTCGAG
# ---------------------- 22222222 --------------- --------------- 11111111 ---------------------------- UUUUUUUUUU TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT


Enh_linker1 = b'GTGA'
Enh_linker2 = b'GACA'

Enh_5p_primer = b"ACAGGAAACTCATGGTGCGT"

Enh_5p_linker1 = b"AATG"
Enh_5p_linker2 = b"CCAC"

Enh_inserts = [b"", b"A", b"GT", b"TCA"]

Tso_capture_seq_Enh_EnhV2 = b"TATGCGTAGTAGGTATG"
Tso_capture_seq_EnhV3 = b"GTGGAGTCGTGATTATA"

#B384_cell_key1 = (
B384_cell_key = [
                [
                  "TGTGTTCGC","TGTGGCGCC","TGTCTAGCG","TGGTTGTCC","TGGTTCCTC","TGGTGTGCT","TGGCGACCG","TGCTGTGGC","TGCTGGCAC","TGCTCTTCC",
                  "TGCCTCACC","TGCCATTAT","TGATGTCTC","TGATGGCCT","TGATGCTTG","TGAAGGACC","TCTGTCTCC","TCTGATTAT","TCTGAGGTT","TCTCGTTCT",
                  "TCTCATCCG","TCCTGGATT","TCAGCATTC","TCACGCCTT","TATGTGCAC","TATGCGGCC","TATGACGAG","TATCTCGTG","TATATGACC","TAGGCTGTG",
                  "TACTGCGTT","TACGTGTCC","TAATCACAT","GTTGTGTTG","GTTGTGGCT","GTTGTCTGT","GTTGTCGAG","GTTGTCCTC","GTTGTATCC","GTTGGTTCT",
                  "GTTGGCGTT","GTTGGAGCG","GTTGCTGCC","GTTGCGCAT","GTTGCAGGT","GTTGCACTG","GTTGATGAT","GTTGATACG","GTTGAAGTC","GTTCTGTGC",
                  "GTTCTCTCG","GTTCTATAT","GTTCGTATG","GTTCGGCCT","GTTCGCGGC","GTTCGATTC","GTTCCGGTT","GTTCCGACG","GTTCACGCT","GTTATCACC",
                  "GTTAGTCCG","GTTAGGTGT","GTTAGAGAC","GTTAGACTT","GTTACCTCT","GTTAATTCC","GTTAAGCGC","GTGTTGCTT","GTGTTCGGT","GTGTTCCAG",
                  "GTGTTCATC","GTGTCACAC","GTGTCAAGT","GTGTACTGC","GTGGTTAGT","GTGGTACCG","GTGGCGATC","GTGCTTCTG","GTGCGTTCC","GTGCGGTAT",
                  "GTGCGCCTT","GTGCGAACT","GTGCAGCCG","GTGCAATTG","GTGCAAGGC","GTCTTGCGC","GTCTGGCCG","GTCTGAGGC","GTCTCAGAT","GTCTCAACC",
                  "GTCTATCGT","GTCGGTGTG","GTCGGAATC","GTCGCTCCG","GTCCTCGCC","GTCCTACCT","GTCCGCTTG","GTCCATTCT","GTCCAATAC","GTCATGTAT",

                  "GTCAGTGGT","GTCAGATAG","GTATTAACT","GTATCAGTC","GTATAGCCT","GTATACTTG","GTATAAGGT","GTAGCATCG","GTACCGTCC","GTACACCTC",
                  "GTAAGTGCC","GTAACAGAG","GGTTGTGTC","GGTTGGCTG","GGTTGACGC","GGTTCGTCG","GGTTCAGTT","GGTTATATT","GGTTAATAC","GGTGTACGT",
                  "GGTGCCGCT","GGTGCATGC","GGTCGTTGC","GGTCGAGGT","GGTAGGCAC","GGTAGCTTG","GGTACATAG","GGTAATCTG","GGCTTGGCC","GGCTTCACG",
                  "GGCTTATGT","GGCTTACTC","GGCTGTCTT","GGCTCTGTG","GGCTCCGGT","GGCTCACCT","GGCGTTGAG","GGCGTGTAC","GGCGTGCTG","GGCGTATCG",
                  "GGCGCTCGT","GGCGCTACC","GGCGAGCCT","GGCGAGATC","GGCGACTTG","GGCCTCTTC","GGCCTACAG","GGCCAGCGC","GGCCAACTT","GGCATTCCT",
                  "GGCATCCGC","GGCATAACC","GGCAACGAT","GGATGTCCG","GGATGAGAG","GGATCTGGC","GGATCCATG","GGATAGGTT","GGAGTCGTG","GGAGAAGGC",
                  "GGACTCCTT","GGACTAGTC","GGACCGTTG","GGAATTAGT","GGAATCTCT","GGAATCGAC","GGAAGCCTC","GCTTGTAGC","GCTTGACCG","GCTTCGGAC",
                  "GCTTCACAT","GCTTAGTCT","GCTGGATAT","GCTGGAACC","GCTGCGATG","GCTGATCAG","GCTGAGCGT","GCTCTTGTC","GCTCTCCTG","GCTCGGTCC",
                  "GCTCCAATT","GCTATTCGC","GCTATGAGT","GCTAGTGTT","GCTAGGATC","GCTAGCACT","GCTACGTAT","GCTAACCTT","GCGTTCCGC","GCGTGTGCC",
                  "GCGTGCATT","GCGTCGGTT","GCGTATGTG","GCGTATACT","GCGGTTCAC","GCGGTCTTG","GCGGCGTCG","GCGGCACCT","GCGCTGGAC","GCGCTCTCC",

                  "GCGCGGCAG","GCGCGATAC","GCGCCGACC","GCGAGCGAG","GCGAGAGGT","GCGAATTAC","GCCTTGCAT","GCCTGCGCT","GCCTAACTG","GCCGTCCGT",
                  "GCCGCTGTC","GCCATGCCG","GCCAGCTAT","GCCAACCAG","GCATGGTTG","GCATCGACG","GCAGGCTAG","GCAGGACGC","GCAGCCATC","GCAGATACC",
                  "GCAGACGTT","GCACTATGT","GCACACGAG","GATTGTCAT","GATTGGTAG","GATTGCACC","GATTCTACT","GATTCGCTT","GATTAGGCC","GATTACGGT",
                  "GATGTTGGC","GATGTTATG","GATGGCCAG","GATCGTTCG","GATCGGAGC","GATCGCCTC","GATCCTCTG","GATCCAGCG","GATACACGC","GAGTTACCT",
                  "GAGTCGTAT","GAGTCGCCG","GAGGTGTAG","GAGGCATTG","GAGCGGACG","GAGCCTGAG","GAGATCTGT","GAGATAATT","GAGACGGCT","GACTTCGTG",
                  "GACTGTTCT","GACTCTTAG","GACCGCATT","GAATTGAGC","GAATATTGC","GAAGGCTCT","GAAGAGACT","GAACTGCCG","GAACGCGTG","CTTGTGTAT",
                  "CTTGTGCGC","CTTGTCATG","CTTGGTCTT","CTTGGTACC","CTTGGATGT","CTTGCTCAC","CTTGCAATC","CTTGAGGCC","CTTGACGGT","CTTCTGATC",
                  "CTTCTCGTT","CTTCTAGGC","CTTCGTTAG","CTTATGTCC","CTTATGCTT","CTTATATAG","CTTAGGTTG","CTTAGGAGC","CTTACTTAT","CTGTTCTCG",
                  "CTGTGCCTC","CTGTCGCAT","CTGTCGAGC","CTGTAGCTG","CTGTACGTT","CTGCTTGCC","CTGCGTAGT","CTGCACACC","CTGATGGAT","CTGAGTCAT",
                  "CTGACGCCG","CTGAACGAG","CTCTTGTAG","CTCTTAGTT","CTCTTACCG","CTCTGCACC","CTCTCGTCC","CTCGTATTG","CTCGACTAT","CTCCTGACG",

                  "CTCACTAGC","CTATACGGC","CGTTCGCTC","CGTTCACCG","CGTATAGTT","CGGTGTTCC","CGGTGTCAG","CGGTCCTGC","CGGCGACTC","CGGCACGGT",
                  "CGGATAGCC","CGGAGAGAT","CGCTAATAG","CGCGTTGGC","CGCGCAGAG","CGCACTGCC","CCTTGTCTC","CCTTGGCGT","CCTTCTGAG","CCTTCTCCT",
                  "CCTTCGACC","CCTTACTTG","CCTGTTCGT","CCTGTATGC","CCTCGGCCG","CCGTTAATT","CCATGTGCG","CCAGTGGTT","CCAGGCATT","CCAGGATCC",
                  "CCAGCGTTG","CATTCCGAT","CATTATACC","CATGTTGAG","ATTGCGTGT","ATTGCGGAC","ATTGCGCCG","ATTGACTTG","ATTCGGCTG","ATTCGCGAG",
                  "ATTCCAAGT","ATTATCTTC","ATTACTGTT","ATTACACTC","ATGTTCTAT","ATGTTACGC","ATGTGTATC","ATGTGGCAG","ATGTCTGTG","ATGGTGCAT",
                  "ATGCTTACT","ATGCTGTCC","ATGCTCGGC","ATGAGGTTC","ATGAGAGTG","ATCTTGGCT","ATCTGTGCG","ATCGGTTCC","ATCATGCTC","ATCATCACT",
                  "ATATCTTAT","ATAGGCGCC","AGTTGGTAT","AGTTGAGCC","AGTGCGACC","AGGTGCTAC","AGGCTTGCG","AGGCCTTCC","AGGCACCTT","AGGAATATG",
                  "AGCGGCCAG","AGCCTGGTC","AGCCTGACT","AGCAATCCG","AGAGATGTT","AGAGAATTC","ACTCGCTTG","ACTCGACCT","ACGTACACC","ACGGATGGT",
                  "ACCAGTCTG","ACATTCGGC","ACATGAGGT","ACACTAATT"
                  ],
                  [
                  "TTGTGTTGT","TTGTGGTAG","TTGTGCGGA","TTGTCTGTT","TTGTCTAAG","TTGTCATAT","TTGTCACGA","TTGTATGAA","TTGTACAGT","TTGGTTAAT",
                  "TTGGTGCAA","TTGGTCGAG","TTGGTATTA","TTGGCACAG","TTGGATACA","TTGGAAGTG","TTGCGGTTA","TTGCCATTG","TTGCACGCG","TTGCAAGGT",
                  "TTGATGTAT","TTGATAATT","TTGAGACGT","TTGACTACT","TTGACCGAA","TTCTGGTCT","TTCTGCACA","TTCTCCTTA","TTCTCCGCT","TTCTAGGTA",
                  "TTCTAATCG","TTCGTCGTA","TTCGTAGAT","TTCGGCTTG","TTCGGAATA","TTCGCCAGA","TTCGATTGT","TTCGATCAG","TTCCTCGGT","TTCCGGCAG",
                  "TTCCGCATT","TTCCAATTA","TTCATTGAA","TTCATGCTG","TTCAGGAGT","TTCACTATA","TTCAACTCT","TTCAACGTG","TTATGCGTT","TTATGATTG",
                  "TTATCCTGT","TTATCCGAG","TTATATTAT","TTAGGCGCG","TTACTGGAA","TTACTAGTT","TTACGTGGT","TTACGATAT","TTACCTAGA","TTACATGAG",
                  "TTACAGCGT","TTACACGGA","TTACACACT","TTAATCAGT","TTAATAGGA","TTAAGTGTG","TTAACCTTG","TTAACACAA","TGTTCACTT","TGTTCAAGA",
                  "TGTTAAGTG","TGTGTTATG","TGTGTCCAA","TGTGGAGCG","TGTCAGTTA","TGTCAGAAG","TGGTTAGTT","TGGTTACAA","TGGCGTTAT","TGGCGCCAA",
                  "TGGAGTCTT","TGCGTATTG","TGATAGAGA","TGAGGTATT","TGAGAATCT","TCTTGGTAA","TCTTCATAG","TCTGTCCTT","TCTGGAATT","TCTACCGCG",
                  "TCGTTCGAA","TCGTCAGTG","TCGACGAGA","TCATGGCTT","TCACACTTA","TATTCCGAA","TATTATGGT","TATGCTATT","TATCAAGGA","TAGTTCAAT",

                  "TAGCTGCTT","TAGAGGAAG","TACCTGTTA","TACACCTGT","GTTGTGCGT","GTTGGCTAT","GTTGCCAAG","GTTGACCTT","GTTCTGCTA","GTTCTGAAT",
                  "GTTCTATCA","GTTCGCGTG","GTTCCTTAT","GTTAGCAGT","GTTACTGTG","GTTACTCAA","GTTAAGAGA","GTTAACTTA","GTGTCGGCA","GTGTCCATT",
                  "GTGCTTGAG","GTGCTCGTT","GTGCTCACA","GTGCCTGGA","GTCTTGTCG","GTCTTGATT","GTCTTCCGT","GTCTTAAGA","GTCTCATCT","GTCTACGAG",
                  "GTCGTTGCT","GTCGTGTTA","GTCGGTAAT","GTCGGATGT","GTCGAGCTG","GTCCGGACT","GTCCAACAT","GTCAGACGA","GTCAGAATT","GTCACTCTT",
                  "GTCAAGGAA","GTATGTCTT","GTATGTACA","GTATCGGTT","GTATATGTA","GTATACAAT","GTAGTTAAG","GTAGTCGAT","GTAGCCTTA","GTAGATACT",
                  "GTACGATTA","GTACAGTCT","GTAATTCGT","GCTTGGCAG","GCTTGCTTG","GCTTGAGGA","GCTTCATTA","GCTTATGCG","GCTGTGTAG","GCTGTCATG",
                  "GCTGGTTGT","GCTGGACTG","GCTGCCTAA","GCTGATATT","GCTCTTAGT","GCTCTATTG","GCTCGCCGT","GCTCCGCTG","GCTATTCTG","GCTATACGA",
                  "GCTACTAAG","GCTACATGT","GCTAACTCT","GCGTTGTAA","GCGTTCTCT","GCGTGCGTA","GCGTCTTGA","GCGTCCGAT","GCGTAAGAG","GCGCTTACG",
                  "GCGCGGATT","GCGCCATAT","GCGCATGAA","GCGATCAAT","GCGAGCCTT","GCGAGATTG","GCGAGAACA","GCCTTGGTA","GCCTTCTAG","GCCTTCACA",
                  "GCCTGAGTG","GCCTCACGT","GCCGGCGAA","GCCGCACAA","GCCATGCTT","GCCATATAT","GCCAATTCG","GCATTCGTT","GCATGATGT","GCAGTTGGA",

                  "GCAGTGTCT","GCACTTGTG","GCAATCTGT","GCAACACTT","GATTGTATT","GATTGCGAG","GATTCCAGT","GATTCATAT","GATTATCAG","GATTAGGTT",
                  "GATGTTGCG","GATGGATCT","GATGCTGAT","GATGCCTTG","GATCTCCTT","GATCGCTTA","GATATTGAA","GATATTACT","GAGTGTTAT","GAGCTCAGT",
                  "GAGCGTGCT","GAGCGTCGA","GAGCGGTTG","GAGCGACTT","GAGCCGAAT","GAGATAGAT","GAGACCTAT","GACGGTCGT","GACGCAGGT","GACGATATG",
                  "GACCTATCT","GAATTAGGA","GAATCAGCT","GAAGTTCAT","GAAGTGGTT","GAAGTATTG","GAAGGCATT","GAACGCTGT","CTTGTCCAG","CTTGGATTG",
                  "CTTGCTGAA","CTTGCCGTG","CTTGATTCT","CTTCTGTCG","CTTCGGCGT","CTTATGAGT","CTTACCGAT","CTGTTAGGT","CTGTCGTCT","CTGTATAAT",
                  "CTGGCTCAT","CTGGATGCG","CTGCGTGTG","CTGCGCGGT","CTGCCGATT","CTGCATTGT","CTGATTAAG","CTGAGATAT","CTGACCTGT","CTCGTATCT",
                  "CTCGGCAAG","CTCGCAATT","CTCCTGCTT","CTCCTAAGT","CTCCGGATG","CTCCGAGCG","CTCACAGGT","CTATTCTAT","CTATTAGTG","CTATGAATT",
                  "CTACATATT","CGTGGCATT","CGTCTTAAT","CGTCTGGTT","CGTCACTGT","CGTAGGTCT","CGGTTCGAG","CGGTTCATT","CGGTGCTCT","CGGTAATTG",
                  "CGGCCTGAT","CGGATATAG","CGGAATATT","CGCTCCAAT","CGCGTTCGT","CGCAGGTTG","CGAGGATGT","CGAGCTGTT","CGACGGCTT","CCTTGTGTG",
                  "CCTGTCTCA","CCTGACTAT","CCTACCTTG","CCGTAGATT","CCGGCTGGT","CATCGGACG","CATCGATAA","CATCCTTCT","CAGTTCTGT","CAGTGCCAG",

                  "CAGGCACTG","CAGCCTCTT","CACTTATAT","CACTGGTCG","CACTGCATG","CACGCGTTG","CACGATGTT","CACCATCTG","CACAGGCGT","ATTGTACAA",
                  "ATTGGTATG","ATTGCTAAT","ATTGCATAG","ATTGCAGTT","ATTCTGCAG","ATTCTACGT","ATTCGGATT","ATTCCGTTG","ATTCATCAA","ATTCAAGAG",
                  "ATTAGCCTT","ATTAATATT","ATGTTAGAG","ATGTTAACT","ATGTAGTCG","ATGGTGTAG","ATGGATTAT","ATCTTGAAG","ATCTGATAT","ATCTCAGAA",
                  "ATCGCTCAA","ATCGCGTCG","ATCCATGGT","ATCATGAGA","ATCATAGTT","ATCAGCGAG","ATCACCATT","ATAGTAATT","ATAGCTGTG","ATACTCTCG",
                  "ATACCTCAT","AGTTGCGCG","AGTTGAATT","AGTTATGAT","AGTGTCCGT","AGTGGCTTG","AGTGCTTCT","AGTATCATT","AGTACACAA","AGGTATGCG",
                  "AGGTATAGT","AGGCTACTT","AGGCCAGGT","AGGAGCGAT","AGCTTATAG","AGCTCTAGA","AGCGTGTAT","AGCGTCACA","AGCCTTCAT","AGCCTGTCG",
                  "AGCCTCGAG","AGCACTGAA","AGATGTACG","AGAGTTAAT","AGACCTCTG","ACTTCTATA","ACTGTCGAG","ACTGTATGT","ACTCTGTAA","ACTCGCGAA",
                  "ACTAGATCT","ACTAACGTT","ACGTTACTG","ACGTGGAAT","ACGGACTCT","ACGCCTAAT","ACGCCGTTA","ACGACGTGT","ACCTCGCAT","ACCATCATA",
                  "ACATATATT","ACAGGCACA","ACACCTGAG","ACACATTCT"
                  ],
                  [
                  "TTGTGGCTG","TTGTGGAGT","TTGTGCGAC","TTGTCTTCA","TTGTAAGAT","TTGGTTCTG","TTGGTGCGT","TTGGTCTAC","TTGGTAACT","TTGGCGTGC",
                  "TTGGATTAG","TTGGAGACG","TTGGAATCA","TTGCGGCGA","TTGCGCTCG","TTGCCTTAC","TTGCCGGAT","TTGCATGCT","TTGCACGTC","TTGCACCAT",
                  "TTGAACCTG","TTCTCGCGT","TTCTCAACT","TTCTACTCA","TTCGTCCAT","TTCGGATAC","TTCGGACGT","TTCGCAATC","TTCCGGTGC","TTCCGACTG",
                  "TTCATTATG","TTCATGGAT","TTCAGCGCA","TTCACCTCG","TTCAAGCAG","TTCAACTAC","TTATGCCAG","TTATGCATC","TTATCGTAC","TTATACCTA",
                  "TTATAATAG","TTATAAGTC","TTAGTTAGC","TTAGCTCAT","TTAGCACTA","TTAGATATG","TTACTACGA","TTACCGTCA","TTACAGAGC","TTAATTGCA",
                  "TTAACAGAT","TGTTGGCTA","TGTTGATGA","TGTTAAGCT","TGTGGCCGA","TGTGCTAGC","TGTGCGTCA","TGTCGCAGT","TGTCGAGCA","TGTACAACG",
                  "TGGTTCCGA","TGGTTCACT","TGGTCAAGT","TGGCTTGTA","TGGCTGTCG","TGGCGTATG","TGGCGCGCT","TGGATGTAC","TGGACTTGC","TGGAATACT",
                  "TGCTAGCGA","TGCGTTGCT","TGCGGTCTG","TGCGCTTAG","TGCGCGACG","TGCCTGCAT","TGCCTAGAC","TGCACGAGT","TGAGTGTGC","TGAGGCTCG",
                  "TCTTCCGTC","TCTTATAGT","TCTTACCAT","TCTGTTGTC","TCTGTTACT","TCTGGCTAG","TCTCAGATC","TCTAGTTGA","TCTAGTACG","TCGTACTAC",
                  "TCGGTGTAG","TCGGCTGCT","TCGCTACTG","TCGATCACG","TCGAGGCAT","TCCGGCGTC","TCCGGAGCT","TCCGCTCGT","TCCGAGTAC","TCCATTCAT",

                  "TCCATGGTC","TCCAAGTCG","TCATTACGT","TCATGCACT","TCAGGTTGC","TCAGACCGT","TCACTCAGT","TCAAGCTCA","TATTGCGCA","TATTCGGCT",
                  "TATTCCAGC","TATTCATCA","TATGTTCAG","TATGGTATG","TATGCAAGT","TATCTGGTC","TATCTGACT","TATCCAGAT","TATCAGTCG","TATCACGCT",
                  "TAGGCGCGA","TAGGCACAT","TAGGATCGT","TAGCATTGC","TAGAGTTAC","TAGACTGAT","TACTTGTCG","TACGTCCGA","TACCGTACT","TACCGCGAT",
                  "TACCAGGAC","TACAGAAGT","TAAGTGCAT","TAAGCTACT","GTTGACCGA","GTTCTCGAC","GTTCCTGCT","GTTATGATG","GTGCTTGCA","GTGCCGCGT",
                  "GTATTGCTG","GTATTCCGA","GTATTAAGC","GTATGACGT","GTAGTTGTC","GTAGTACAT","GTAGCTCGA","GGTTGCTCA","GGTTGAGTA","GGTTAACGT",
                  "GGTGTGGCA","GGTCTTCAG","GGTCGTCTA","GGTCGGCGT","GGTCCGACT","GGTCATGTC","GGTCACATG","GGTAGTGCT","GGTAGCGTC","GGTACCAGT",
                  "GGTAAGGAT","GGCTTGTGC","GGCTTGACT","GGCTTACGA","GGCTGTAGT","GGCTGGCAG","GGCTCCATC","GGCGTGGAT","GGCGTAATC","GGCGCAAGT",
                  "GGCGAGTAG","GGCGACCGT","GGCCTGTCA","GGCCATTGC","GGCACTCTG","GGATGTCAT","GGAGTAACT","GGAGAACGA","GGACTGGCT","GGACGTTCA",
                  "GGAACGTGC","GCTGTCCAT","GCTGGTTCA","GCTGCAACT","GCTCGTTAC","GCTATAGAT","GCTAGTCGT","GCTACCATG","GCGTTCTGA","GCGTGTTAG",
                  "GCGGTATCG","GCGGAGCAT","GCGCGGTGC","GCGCCTAGT","GCGCCGGCT","GCCTTCATG","GCCATACTG","GCATGTTGA","GCATGCTAC","GCAGTATAC",

                  "GCAGGTACT","GCAGCGCGT","GCACCTCAT","GCAATTCGA","GATTGCCGT","GATGAACAT","GATCTTCGA","GATCTGCAT","GAGTGGCAT","GAGTCGGAC",
                  "GAGTATGAT","GAGGCGAGT","GAGGCAACG","GAGCGCACT","GAATAGGCT","ATTGTCACT","ATTGTATCA","ATTGGTCAG","ATTGGCGAT","ATTGATCGT",
                  "ATTCGTAGT","ATTCATACG","ATTCAGGAC","ATTACTTCA","ATTAATTAG","ATTAAGCAT","ATGTCTCTA","ATGTAGCGT","ATGGCATAC","ATGGAGATC",
                  "ATGGACTCG","ATGGAACGA","ATGCTTCAT","ATGCTCGCT","ATGCGACGT","ATGCCGTAG","ATGAGTTCG","ATGACTATC","ATGACCGAC","ATCTTATGC",
                  "ATCTTACTA","ATCTATCAG","ATCGTGTAC","ATCGTCTGA","ATCGGCATG","ATCGCGAGC","ATCGCAACG","ATCGATGCT","ATCGAATAG","ATCCTTCTG",
                  "ATCCTGCGT","ATCCGCACT","ATCCATTAC","ATCCAAGCA","ATCAGATCA","ATCACACAT","ATCAACGTC","ATCAACCGA","ATATTGAGT","ATATTCGTC",
                  "ATATTACAG","ATATCTTGA","ATATCGCAT","ATATCAATC","ATAGTCCTG","ATAGGTCTA","ATAGCTGAC","ATAGCGGTA","AGTTCGCTG","AGTTACAGC",
                  "AGTTAACTA","AGTGCAATC","AGTCTGGTA","AGTCTGAGC","AGTCTACAT","AGTCGAACT","AGTCCATCG","AGTCATTCA","AGTATCCAG","AGTAGACTG",
                  "AGTAATCGA","AGTAAGTGC","AGGTTGGCT","AGGTTCTAG","AGGTGTTCA","AGGTGCCAT","AGGTCTGAT","AGGTCGTAC","AGGTCAGCA","AGGCTTATC",
                  "AGGCTATGA","AGGCCGACG","AGGCCAAGC","AGGCAGGTC","AGGCAAGAT","AGGAGCAGT","AGGACCGCT","AGGAATTAC","AGCTTGGAC","AGCTTAAGT",

                  "AGCTACACG","AGCGTTACG","AGCGGTGCA","AGCGGAGTC","AGCGGACGA","AGCGCGCTA","AGCGATAGC","AGCGACTCA","AGCCTCTAC","AGCCGTCGT",
                  "AGCATGATC","AGCACTTCG","AGCACGGCA","AGATTCTGA","AGATTAGAT","AGATGATAG","AGATATGTA","AGATACCGT","AGAGTGCGT","AGAGCCGAT",
                  "AGACTCACT","ACTTGCCTA","ACTTGAGCA","ACTTCTAGC","ACTTCGACT","ACTTAGTAC","ACTGTTGAT","ACTGTAACG","ACTGGTATC","ACTGACGTC",
                  "ACTGAAGCT","ACTCTGATG","ACTCCTGAC","ACTCCGCTA","ACTCAACTG","ACTATTGCA","ACTAGGCAG","ACTACGCGT","ACTAATACT","ACGTTCGTA",
                  "ACGTGTGCT","ACGTGTATG","ACGTGGAGC","ACGTCTTCG","ACGTCAGTC","ACGGTCTCA","ACGGTCCGT","ACGGTACAG","ACGGCGCTG","ACGCTGCGA",
                  "ACGCGTGTA","ACGCGCCAG","ACGATGTCG","ACGATGGAT","ACGATCTAC","ACGAGCTGA","ACGAGCATC","ACGAATCGT","ACGAACGCA","ACCTTGTAG",
                  "ACCTGTTGC","ACCTGTCAT","ACCTCGATC","ACCTAGGTA","ACCTACTGA","ACCTAATCG","ACCGTAGCA","ACCGGTAGT","ACCGGCTAC","ACCGCTTCA",
                  "ACATTGTGC","ACATTCTCG","ACATGGCTG","ACATGACGA","ACATATGAT","ACATATACG","ACAGCGTAC","ACACTTGCT","ACACTATCA","ACACGCATG",
                  "ACACCAGTA","ACACCAACT","ACACATAGT","ACACACCTA"
                  ]
                  ]

for x in range(2):
    #convert to binary
    B384_cell_key[x] = [x.encode('ascii') for x in B384_cell_key[x]]
    

def hamming(a, b):
    la = len(a)
    lb = len(b)
    if la != lb:
        return la
    return sum([a[x] != b[x] for x in range(len(a))])

def correct_bc(bc, bc_list, distance=1):
    if bc in bc_list:
        return bc
    distances = np.array([hamming(bc, sh_bc) for sh_bc in bc_list])
    accepted = np.where(distances <= distance)[0]
    if len(accepted) >0:
        # return the first one
        return bc_list[accepted[0]]
    else:
        return b''
 

def get_options():
    parser = argparse.ArgumentParser(prog='process_share.py')
    parser.add_argument('-i', '--cls_read', help='Fastq file with cell barcodes. R1 for scRNA, I1 for scATAC')
    parser.add_argument('-R', '--rna', help='Process as scRNA-seq', action='store_true')
    parser.add_argument('-L', '--umi_length', help='Length of UMI (stiched to CB)', default=8)
    parser.add_argument('-p', '--prefix', help='Prefix for output files')
    parser.add_argument('-C', '--correct_barcodes', help='Correct barcode to expected ones', action='store_true')
    parser.add_argument('-t', '--threshold', help='Max distance when correcting barcodes', default=1, type=int)
    parser.add_argument('-n', '--n_seq', help='Max number of sequences to process (for debugging)', default=0, type=int)
  
    options = parser.parse_args()
  
    return options

def main():
    nl = b'\n'
    dnl = b'\n+\n'
    dark = b'GGGGGGGGGGGGGGGGGGGG'
    _chunk_size = 512 # number 

    options = get_options()
    
    options.correct_barcodes 
    
    if options.rna:
        spacer_1, spacer_2 = Enh_linker1, Enh_linker2
    else:
        spacer_1, spacer_2 = Enh_5p_linker1, Enh_5p_linker2

    read_iterator = HTSeq.FastqReader(options.cls_read)
    
    fout1 = f'{options.prefix}_RB.fastq.gz' # for Barcode
    
    raw_out1 = open(fout1, 'wb')
    fh_out1 = bgzip.BGZipWriter(raw_out1, batch_size=256)

    r1_spool = b''
    _spool_counter = 0
    
    n_tot = 0
    n_spwrong = 0
    n_fail = 0
    
    for item in tqdm(read_iterator):
    
        if options.n_seq > 0 and n_tot == options.n_seq:
            break
    
        seq1 = item.seq
        qual1 = item.qualstr

        cls1 = cls2 = cls3 = b'NNNNNNNNN'
        q_cls1 = q_cls2 = q_cls3 = b'!!!!!!!!!'
        umi = b'NNNNNNNN'
        q_umi = b'!!!!!!!!'
        
        if options.rna:
            # isolate 3 barcodes and UMI
            # spacer_1 and spacer_2 are at variable position, depending on
            # 5' diversity
            it_c = 0
            for x in range(4):
                if seq1.startswith(Enh_inserts[x]) and seq1[x+9:x+13] == spacer_1 and seq1[x+22:x+26] == spacer_2:
                    cls1, cls2, cls3, umi = (seq1[x:x+9], seq1[x+13:x+22], seq1[x+26:x+35], seq1[x+35:x+43] )
                    q_cls1, q_cls2, q_cls3, q_umi = (qual1[x:x+9], qual1[x+13:x+22], qual1[x+26:x+35], qual1[x+35:x+43] )
                    it_c += 1
            if it_c == 0:
                n_spwrong += 1

        else:
            # isolate 3 barcodes, removing 5p linker
            # spacer_1, spacer_2 are at fixed
            if hamming(seq1[:20], Enh_5p_primer) < options.threshold and seq1[29:33] == spacer_1 and seq1[42:46] == spacer_2:
                cls1, cls2, cls3 = seq1[20:29], seq1[33:42], seq1[46:55]
                q_cls1, q_cls2, q_cls3 = qual1[20:29], qual1[33:42], qual1[46:55]
            else:
                n_spwrong += 1


        n_tot += 1

        name1 = bytes('@' + item.name, encoding='ascii')

        if options.correct_barcodes:
            cls1 = correct_bc(cls1, B384_cell_key[0], options.threshold)
            cls2 = correct_bc(cls2, B384_cell_key[1], options.threshold)
            cls3 = correct_bc(cls3, B384_cell_key[2], options.threshold)
            
        seq1_out = cls1 + cls2 + cls3
        q_seq1_out = q_cls1 + q_cls2  + q_cls3
        
        if options.rna:
            seq1_out = seq1_out + umi
            q_seq1_out = q_seq1_out + q_umi
        
        # since bccorrection returns empty bc if not found
        # we can use it to skip bad reads
        
        if options.correct_barcodes  and len(seq1_out) != len(q_seq1_out):
            n_fail += 1
            seq1_out = b'NNNNNNNNNNNNNNNNNNNNNNNNNNN'
        
        r1_spool = r1_spool + name1 + nl + seq1_out + dnl + q_seq1_out + nl
        
        _spool_counter += 1
        
        if _spool_counter == _chunk_size:
            fh_out1.write(r1_spool)
            r1_spool = b''
            _spool_counter = 0
    
    # end, write remaining spool and close files
    fh_out1.write(r1_spool)
    fh_out1.close()
    raw_out1.close()

    n_pass = n_tot - n_spwrong - n_fail
    sys.stderr.write(f"Total sequences:\t{n_tot}\n")
    f = n_spwrong / n_tot * 100
    sys.stderr.write(f"Error in spacers:\t{n_spwrong} ({f:.3f}%)\n")
    f = n_fail / (n_tot - n_spwrong) * 100
    sys.stderr.write(f"Failed BC:\t{n_fail} ({f:.3f}%)\n")
    f = n_pass / n_tot * 100
    sys.stderr.write(f"Passing sequences:\t{n_pass} ({f:.3f}%)\n")
    
#    eff = n_pass / n_tot * 100
#    sys.stderr.write(f'Found {n_pass} out of {n_tot} sequences {eff:.3f}%\n')
#    eff = n_spwrong / n_tot * 100
#    sys.stderr.write(f'Found {n_spwrong} dark sequences {eff:.3f}%\n')
#    eff = n_fail / (n_tot - n_spwrong) * 100
#    sys.stderr.write(f'Could not fix barcode for {n_fail} sequences {eff:.3f}%\n')



if __name__ == '__main__':
    main()

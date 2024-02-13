#!/usr/bin/env python3

##########################################################################################################
#### Merge compiled TCGA mRNA fpkm table with TCGA isomiRNA RPM table 
#### JToubia - January 2021
#########

### import modules
import sys
import numpy as np
import pandas as pd

### Read in arguments
rnaseq_summary = sys.argv[1]
mirseq_summary = sys.argv[2]
tcga_ID = sys.argv[3]

### Merge the dfs and save
rnaseq_df = pd.read_csv(rnaseq_summary, sep='\t', index_col=0, header=0)
mirseq_df = pd.read_csv(mirseq_summary, sep='\t', index_col=0, header=0)

rnaseq_df = rnaseq_df.transpose()
mirseq_df = mirseq_df.transpose()
mirseq_df.drop(mirseq_df.index[:1], inplace=True)

final_df = rnaseq_df.join(mirseq_df, how='outer')
final_df = final_df[~final_df.index.duplicated(keep='first')]
final_df.to_csv('GDC_' + tcga_ID + '_gene_mir_exp.tsv', sep='\t')

final_df = final_df.replace(0, np.nan)
final_df.to_csv('GDC_' + tcga_ID + '_gene_mir_exp_zero2nan.tsv', sep='\t')






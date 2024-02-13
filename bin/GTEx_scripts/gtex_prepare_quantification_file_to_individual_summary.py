#!/usr/bin/env python3

##########################################################################################################
#### Convert GTEx multi patient/tissue type expression file to dataframes and tables
#### JToubia - January 2021
#########

### import modules
import os
import sys
import pandas as pd
import numpy as np
from pathlib import Path

### Read in arguments
sam_attr = sys.argv[1]
exp_datafile = sys.argv[2]
date_stamp = sys.argv[3]
data_type = sys.argv[4]

### Read in sample attributes
sam_attr_data = pd.read_csv(sam_attr , sep='\t', index_col=0)

### Retain only the columns we need
sam_attr_data_cut = sam_attr_data[['SMTS', 'SMAFRZE']]
sam_attr_data_cut = sam_attr_data_cut.loc[sam_attr_data_cut['SMAFRZE'] == "RNASEQ"]

### Get only unique elements for tissue type
sam_attr_data_cut['SMTS'] = sam_attr_data_cut['SMTS'].str.replace(' ', '_')
SMTS_unique_list = list(sam_attr_data_cut['SMTS'].unique())
SMTS_unique_list_len = str(len(SMTS_unique_list))

### Create a directory to hold each tissue types data
for dirname in SMTS_unique_list:
	Path(dirname).mkdir(parents=True, exist_ok=True)

### Subset the data by tissue type and save
group_by_sam_attr_data_cut = sam_attr_data_cut.groupby(sam_attr_data_cut['SMTS'])

index_number = 1
for key, values in group_by_sam_attr_data_cut: 
	print("\tsubsetting data for " + key + ": Dataset " + str(index_number) + " out of " + SMTS_unique_list_len)
	cols_to_get = list(values.index)
	cols_to_get = ["Description"] + cols_to_get
	exp_data = pd.read_csv(exp_datafile , sep='\t', index_col=0, skiprows=2, usecols=cols_to_get)
	### Do final cleaning
	exp_data = exp_data[~exp_data.index.duplicated(keep='first')]
	exp_data = exp_data.replace(0, np.nan)
	### Save to file
	exp_data.to_csv(key + '/' + "GTEx-" + key + "-" + date_stamp + "_" + data_type + '-mRNA.tsv', sep='\t')
	index_number += 1




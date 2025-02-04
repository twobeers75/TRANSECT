#!/usr/bin/env python3

##########################################################################################################
#### Rename raw TCGA expression files from "file.id" to "patient.id"
#### JToubia - January 2021
#########

### import modules
import sys
import pandas as pd
import os

### Read in arguments
sample_sheet = sys.argv[1]
file_extension_orig = sys.argv[2]
file_extension_new = sys.argv[3]

### Do the renaming
sample_sheet_df = pd.read_csv(sample_sheet, sep="\t")
for index, row_data in sample_sheet_df.iterrows():
	if file_extension_orig[-3:] == ".gz":
		filename = row_data["file_name"][:-3]
	else:
		filename = row_data["file_name"]
		
	SampleID_filename = row_data["cases.0.samples.0.submitter_id"] + file_extension_new
	os_mv_command = "mv " + filename + " " + SampleID_filename
	os.system(os_mv_command)

### prepare mock manifest file
outfile = os.path.join(os.path.dirname(sample_sheet), "MOCK_MANIFEST.txt")
mock_manifest = sample_sheet_df[["id","file_name"]]
mock_manifest.to_csv(outfile, sep='\t', index=False)


#!/usr/bin/env python3

##########################################################################################################
#### Download TCGA sample sheet for custom analyses
#### USAGE example: Download_SampleSheet.py TCGA-OV mRNA_expression_counts
#### JToubia - January 2021
#########

import sys
import re
import requests
import json

### for RNA-seq data use: field: files.analysis.workflow_type with value: [STAR - Counts]
### for miRNA-seq data use: field: files.data_type with value: [miRNA Expression Quantification or Isoform Expression Quantification]

project_id = sys.argv[1] #project_id = "TCGA-BRCA"
folder_name = sys.argv[2] #folder_name = "mRNA_expression_counts"

files_endpt = "https://api.gdc.cancer.gov/files"

if "mRNA" in folder_name:
	field_name = "files.analysis.workflow_type"
	value_name = "STAR - Counts"
	filters = {
		"op": "and",
		"content":[
			{
			"op": "in",
			"content":{
				"field": "cases.project.project_id",
				"value": [project_id]
				}
			},
			{
			"op": "in",
			"content":{
				"field": "files.data_type",
				"value": ["Gene Expression Quantification"]
				}
			},
			{
			"op": "in",
			"content":{
				"field": field_name,
				"value": [value_name]
				}
			}
		]
	}
else:
	field_name = "files.data_type"
	if "isomiR" in folder_name:
		value_name = "Isoform Expression Quantification"
	else: 
		value_name = "miRNA Expression Quantification"
	filters = {
		"op": "and",
		"content":[
			{
			"op": "in",
			"content":{
				"field": "cases.project.project_id",
				"value": [project_id]
				}
			},
			{
			"op": "in",
			"content":{
				"field": field_name,
				"value": [value_name]
				}
			}
		]
	}

# define the fields we need (NOTE: not all these fields are needed but they are useful)
fields = [
	"file_id",
	"file_name",
	"data_type",
	"cases.project.project_id",
	"cases.samples.submitter_id",
	"cases.samples.sample_type"
	]

fields = ",".join(fields)

# A POST is used, so the filter parameters can be passed directly as a Dict object.
params = {
	"filters": json.dumps(filters),
	"fields": fields,
	"format": "TSV",
	"size": "5000"
	}
 
# The parameters are passed to 'json' rather than 'params' in this case
response = requests.post(files_endpt, headers = {"Content-Type": "application/json"}, json = params)

#print(response.content.decode("utf-8"))
with open("sample_sheet.tsv", "wb") as output_file:
	output_file.write(response.content)



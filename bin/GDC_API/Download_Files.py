#!/usr/bin/env python3

##########################################################################################################
#### Download TCGA expression files for custom analyses
#### USAGE example: Download_Files.py TCGA-OV mRNA_expression_counts
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

# Here a GET is used, so the filter parameters should be passed as a JSON string.
params = {
	"filters": json.dumps(filters),
	"fields": "file_id",
	"format": "JSON",
	"size": "10000"
	}

response = requests.get(files_endpt, params = params)

file_uuid_list = []

# This step populates the download list with the file_ids from the previous query
for file_entry in json.loads(response.content.decode("utf-8"))["data"]["hits"]:
	file_uuid_list.append(file_entry["file_id"])

data_endpt = "https://api.gdc.cancer.gov/data"

params = {"ids": file_uuid_list}

response = requests.post(data_endpt, data = json.dumps(params), headers = {"Content-Type": "application/json"})

response_head_cd = response.headers["Content-Disposition"]

file_name = re.findall("filename=(.+)", response_head_cd)[0]

with open(file_name, "wb") as output_file:
	output_file.write(response.content)



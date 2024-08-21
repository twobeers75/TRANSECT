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
import glob
import os
from os import system

### for RNA-seq data use: field: files.analysis.workflow_type with value: [STAR - Counts]
### for miRNA-seq data use: field: files.data_type with value: [miRNA Expression Quantification or Isoform Expression Quantification]

project_id = sys.argv[1] #project_id = "TCGA-BRCA"
folder_name = sys.argv[2] #folder_name = "mRNA_expression_counts" or folder_name = "isomiR_expression_rpm"

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

uuid_list_len = len(file_uuid_list)
print("\t" + str(uuid_list_len) + " uuids for " + project_id + " " + folder_name)

data_endpt = "https://api.gdc.cancer.gov/data"

# Splitting post request (since April 2024) to bypass "Content-Disposition" errors when downloading big datasets (>500 uuids)
for i in range(0, uuid_list_len, 100):
	file_start = i
	file_end = i + 100
	if (file_end > uuid_list_len):
		file_end = uuid_list_len
	
	print("\tretrieving files " + str(file_start+1) + " to " + str(file_end) + " of " + str(uuid_list_len) + " in uuid list")
	params = {"ids": file_uuid_list[file_start:file_end]}
	response = requests.post(data_endpt, data = json.dumps(params), headers = {"Content-Type": "application/json"})
	response_head_cd = response.headers["Content-Disposition"]
	file_name = re.findall("filename=(.+)", response_head_cd)[0]
	
	with open(file_name, "wb") as output_file:
		output_file.write(response.content)

# Now need to combine multiple gz files (if present) before continuing
gz_files = glob.glob("*.tar.gz")
if (len(gz_files) > 1):
	print("\tcombining " + str(len(gz_files)) + " downloaded packages")
	system("cat gdc_download_*.tar.gz > gdc_download_combined.tar.gz")

# Remove original gz files if they were combined
if (len(gz_files) > 1):
	for file in gz_files:
		os.remove(file)







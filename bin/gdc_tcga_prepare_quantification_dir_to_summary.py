#!/usr/bin/env python3

##########################################################################################################
#### Convert TCGA single patient expression files to dataframes and tables
#### JToubia - January 2021
#########

### import modules
import os
import sys
import glob
import pandas as pd
from os.path import basename
from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

### Read in arguments
exp_dir_name = sys.argv[1]
tcga_ID = sys.argv[2]
rna_type = sys.argv[3]
ref_file_folder = sys.argv[4]
non_tcga_data = sys.argv[5]

### Functions
def create_template_dfs(main_col_name, lookup_df_filename):
	### set column names
	column_names = [main_col_name, 'count>0', 'mean', 'std', 'min', 'max']
	### Read in the lookup table
	lookup_df = pd.read_csv(lookup_df_filename , sep='\t', index_col=-0, header=None, names=[main_col_name])
	### create blank dfs
	df1 = pd.DataFrame(data=None, index=None, columns=column_names)
	df2 = pd.DataFrame(data=None, index=None, columns=column_names)
	df3 = pd.DataFrame(data=None, index=None)
	### Add Ensembl ID's as index and gene names
	df1[main_col_name] = lookup_df[main_col_name]
	df2[main_col_name] = lookup_df[main_col_name]
	df3[main_col_name] = lookup_df[main_col_name]
	
	return df1, df2, df3

def get_sample_code(file):
	filename = basename(file)
	TCGA_barcode = filename.split('.')[0]
	project,tss,participant,sample_vial = TCGA_barcode.split("-")
	sample_code = sample_vial[:2]
	
	return TCGA_barcode, sample_code

def get_sample_code_non_TCGA(file):
	filename = basename(file)
	GDC_barcode = filename.split('.')[0]
	
	return GDC_barcode

def get_simple_stats(df):
	dfcount = df.count(axis=1)-1
	dfmean = df.mean(axis=1, numeric_only=True)
	dfstd = df.std(axis=1, numeric_only=True)
	dfmin = df.min(axis=1, numeric_only=True)
	dfmax = df.max(axis=1, numeric_only=True)
	df['count>0'] = dfcount
	df['mean'] = dfmean
	df['std'] = dfstd
	df['min'] = dfmin
	df['max'] = dfmax
	df = df.sort_values('mean', ascending=False, axis=0)
	
	return df

### Create dict of all the TCGA sample types
SAMPLE_TYPE_CANCER = {'01': 'Primary solid Tumor',
	'02': 'Recurrent Solid Tumor',
	'03': 'Primary Blood Derived Cancer - Peripheral Blood',
	'04': 'Recurrent Blood Derived Cancer - Bone Marrow',
	'05': 'Additional - New Primary',
	'06': 'Metastatic',
	'07': 'Additional Metastatic',
	'08': 'Human Tumor Original Cells',
	'09': 'Primary Blood Derived Cancer - Bone Marrow'}
SAMPLE_TYPE_NORMAL = {'10': 'Blood Derived Normal',
	'11': 'Solid Tissue Normal',
	'12': 'Buccal Cell Normal',
	'13': 'EBV Immortalized Normal',
	'14': 'Bone Marrow Normal'}
SAMPLE_TYPE_OTHER = {'20': 'Control Analyte',
	'40': 'Recurrent Blood Derived Cancer - Peripheral Blood',
	'50': 'Cell Lines',
	'60': 'Primary Xenograft Tissue',
	'61': 'Cell Line Derived Xenograft Tissue'}

### Process the data
### Note: these blocks are very similar but enough different that I have not created a single function
if rna_type == 'mRNA':
	### Create template df's that we will need 
	lookup_df_filename = ref_file_folder + '/gencode.vXX.annotation.lookup'
	tumor_df, normal_df, all_df = create_template_dfs('gene_name', lookup_df_filename)
	fpkm_tumor_df, fpkm_normal_df, fpkm_all_df = create_template_dfs('gene_name', lookup_df_filename)
	tpm_tumor_df, tpm_normal_df, tpm_all_df = create_template_dfs('gene_name', lookup_df_filename)
	
	### Start processing each of the individual files
	mrna_files = glob.glob(exp_dir_name) 
	if non_tcga_data == 'true':
		for file in mrna_files:
			GDC_barcode = get_sample_code_non_TCGA(file)
			
			temp_df = pd.read_csv(file, sep='\t', index_col=0, header=1)
			tumor_df[GDC_barcode] = temp_df['unstranded']
			all_df[GDC_barcode] = temp_df['unstranded']
			fpkm_all_df[GDC_barcode] = temp_df['fpkm_unstranded']
			tpm_all_df[GDC_barcode] = temp_df['tpm_unstranded']
	else:
		for file in mrna_files:
			TCGA_barcode, sample_code = get_sample_code(file)
			
			if sample_code in SAMPLE_TYPE_CANCER.keys():
				temp_df = pd.read_csv(file, sep='\t', index_col=0, header=1)
				tumor_df[TCGA_barcode] = temp_df['unstranded']
				all_df[TCGA_barcode] = temp_df['unstranded']
				fpkm_all_df[TCGA_barcode] = temp_df['fpkm_unstranded']
				tpm_all_df[TCGA_barcode] = temp_df['tpm_unstranded']
			elif sample_code in SAMPLE_TYPE_NORMAL.keys():
				temp_df = pd.read_csv(file, sep='\t', index_col=0, header=1)
				normal_df[TCGA_barcode] = temp_df['unstranded']
				all_df[TCGA_barcode] = temp_df['unstranded']
				fpkm_all_df[TCGA_barcode] = temp_df['fpkm_unstranded']
				tpm_all_df[TCGA_barcode] = temp_df['tpm_unstranded']
			elif sample_code in SAMPLE_TYPE_OTHER.keys():
				print('Skipping Other Sample Type with TCGA bacode: ' + TCGA_barcode)
			else:
				print('Unknown samplecode (' + sample_code + ') in TCGA bacode: ' + TCGA_barcode)
	
	### a significant number of entries with duplicate gene_names (mostly smallRNAs and ncRNAs) so need to dedup (keeping the highest expressed).
	all_df['mean'] = all_df.mean(numeric_only=True, axis=1)
	all_df = all_df.sort_values(by=['gene_name', 'mean'])
	all_df = all_df.drop_duplicates(subset='gene_name', keep='last')
	all_df = all_df.drop(['mean'], axis=1)
	
	fpkm_all_df['mean'] = fpkm_all_df.mean(numeric_only=True, axis=1)
	fpkm_all_df = fpkm_all_df.sort_values(by=['gene_name', 'mean'])
	fpkm_all_df = fpkm_all_df.drop_duplicates(subset='gene_name', keep='last')
	fpkm_all_df = fpkm_all_df.drop(['mean'], axis=1)
	
	tpm_all_df['mean'] = tpm_all_df.mean(numeric_only=True, axis=1)
	tpm_all_df = tpm_all_df.sort_values(by=['gene_name', 'mean'])
	tpm_all_df = tpm_all_df.drop_duplicates(subset='gene_name', keep='last')
	tpm_all_df = tpm_all_df.drop(['mean'], axis=1)
	
	### write the normalised outputs here
	if non_tcga_data == 'true':
		fpkm_all_df.to_csv('GDC_' + tcga_ID + '_FPKM-mRNA_all.tsv', sep='\t', index=False)
		tpm_all_df.to_csv('GDC_' + tcga_ID + '_TPM-mRNA_all.tsv', sep='\t', index=False)
	else:
		fpkm_all_df.to_csv('GDC_' + tcga_ID.replace("Count", "FPKM") + '_all.tsv', sep='\t', index=False)
		tpm_all_df.to_csv('GDC_' + tcga_ID.replace("Count", "TPM") + '_all.tsv', sep='\t', index=False)
	
elif rna_type == 'miRNA':
	### Create template df's that we will need 
	lookup_df_filename = ref_file_folder + '/miRname_tcga.lookup'
	tumor_df, normal_df, all_df = create_template_dfs('miR_name', lookup_df_filename)
	
	### Start processing each of the individual files
	mir_files = glob.glob(exp_dir_name)
	for file in mir_files:
		TCGA_barcode, sample_code = get_sample_code(file)
		
		if sample_code in SAMPLE_TYPE_CANCER.keys():
			temp_df = pd.read_csv(file, sep='\t', index_col=0)
			tumor_df[TCGA_barcode] = temp_df['reads_per_million_miRNA_mapped']
			all_df[TCGA_barcode] = temp_df['reads_per_million_miRNA_mapped']
		elif sample_code in SAMPLE_TYPE_NORMAL.keys():
			temp_df = pd.read_csv(file, sep='\t', index_col=0)
			normal_df[TCGA_barcode] = temp_df['reads_per_million_miRNA_mapped']
			all_df[TCGA_barcode] = temp_df['reads_per_million_miRNA_mapped']
		elif sample_code in SAMPLE_TYPE_OTHER.keys():
			print('Skipping Other Sample Type with TCGA bacode: ' + TCGA_barcode)
		else:
			print('Unknown samplecode (' + sample_code + ') in TCGA bacode: ' + TCGA_barcode)
	
elif rna_type == 'isomiRNA': # this one is complicated!
	### Create template df's that we will need 
	lookup_df_filename = ref_file_folder + '/miRBase_21.lookup'
	tumor_df, normal_df, all_df = create_template_dfs('miR_name', lookup_df_filename)
	
	### Start processing each of the individual files
	mir_files = glob.glob(exp_dir_name)
	num_files = len(mir_files)
	for i, file in enumerate(mir_files):
		#for testing
		#if i == 20:
		#	break
		if i % 100 == 0:
			print("processing file: " + str(i) + " of " + str(num_files))
		
		TCGA_barcode, sample_code = get_sample_code(file)
		
		if sample_code in SAMPLE_TYPE_CANCER.keys():
			raw_df = pd.read_csv(file, sep='\t', index_col=0)
			groupby_df = raw_df.groupby(raw_df['miRNA_region'])
			temp_df = pd.DataFrame(data=None, index=None, columns=['MIMAT_ID', 'reads_per_million_miRNA_mapped'])
			temp_list = []
			for key, values in groupby_df: 
				#temp_df = temp_df.append({'MIMAT_ID': key, 'reads_per_million_miRNA_mapped': sum(values['reads_per_million_miRNA_mapped'])}, ignore_index=True)
				temp_dict = {'MIMAT_ID': key, 'reads_per_million_miRNA_mapped': sum(values['reads_per_million_miRNA_mapped'])}
				temp_list.append(temp_dict)
			temp_df = temp_df.from_records(temp_list)
			temp_df = temp_df.set_index('MIMAT_ID')
			temp_df.columns = [TCGA_barcode]
			tumor_df = pd.concat([tumor_df, temp_df], axis=1, sort=False)
			all_df = pd.concat([all_df, temp_df], axis=1, sort=False)
		elif sample_code in SAMPLE_TYPE_NORMAL.keys():
			raw_df = pd.read_csv(file, sep='\t', index_col=0)
			groupby_df = raw_df.groupby(raw_df['miRNA_region'])
			temp_df = pd.DataFrame(data=None, index=None, columns=['MIMAT_ID', 'reads_per_million_miRNA_mapped'])
			temp_list = []
			for key, values in groupby_df: 
				#temp_df = temp_df.append({'MIMAT_ID': key, 'reads_per_million_miRNA_mapped': sum(values['reads_per_million_miRNA_mapped'])}, ignore_index=True)
				temp_dict = {'MIMAT_ID': key, 'reads_per_million_miRNA_mapped': sum(values['reads_per_million_miRNA_mapped'])}
				temp_list.append(temp_dict)
			temp_df = temp_df.from_records(temp_list)
			temp_df = temp_df.set_index('MIMAT_ID')
			temp_df.columns = [TCGA_barcode]
			normal_df = pd.concat([normal_df, temp_df], axis=1, sort=False)
			all_df = pd.concat([all_df, temp_df], axis=1, sort=False)
		elif sample_code in SAMPLE_TYPE_OTHER.keys():
			print('Skipping Other Sample Type with TCGA bacode: ' + TCGA_barcode)
		else:
			print('Unknown samplecode (' + sample_code + ') in TCGA bacode: ' + TCGA_barcode)
	
	### clean up all dfs of unwanted rows
	if "precursor" in tumor_df.index:
		tumor_df = tumor_df.drop(["precursor","unannotated","stemloop"])
	if "precursor" in normal_df.index:
		normal_df = normal_df.drop(["precursor","unannotated","stemloop"])
	if "precursor" in all_df.index:
		all_df = all_df.drop(["precursor","unannotated","stemloop"])
	
else:
	print('Unknown RNA species: ' + rna_type)
	print('Warning: skipping this data') 

### All done, now finish up
### Add simple stats to tables
final_tumor_df = get_simple_stats(tumor_df)
final_normal_df = get_simple_stats(normal_df)

### Write out the final tables
final_tumor_df.to_csv('GDC_' + tcga_ID + '_tumor.tsv', sep='\t')
final_normal_df.to_csv('GDC_' + tcga_ID + '_normal.tsv', sep='\t')
all_df.to_csv('GDC_' + tcga_ID + '_all.tsv', sep='\t', index=False)







#!/usr/bin/env bash

##########################################################################################################
#### Prepare GTEx data for custom analyses
#### JToubia - January 2021
#########

##########################################################################################################
#### Source setup
#########
SCRIPT_FOLDER=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
BASE_FOLDER=$(dirname ${SCRIPT_FOLDER})
DATA_FOLDER="${BASE_FOLDER}/data/GTEx"
REF_FILES_FOLDER="${BASE_FOLDER}/REF_FILES"

##########################################################################################################
#### Functions
#########
process_folder( ) {
	date_stamp=$1
	sam_attr=$2
	exp_data=$3
	data_type=$4
	
	if (${extract_flag}); then extract_files; fi
	create_ind_summary
}

extract_files( ) {
	echo "unziping files"
	cd ${working_folder}
	for file in ./*.gz; do 
		gunzip ${file}
	done
}

create_ind_summary( ) {
	echo "creating individual summaries"
	cd ${working_folder}
	${SCRIPT_FOLDER}/gtex_prepare_quantification_file_to_individual_summary.py \
		"${working_folder}/${sam_attr}" \
		${exp_data} ${date_stamp} ${data_type}
}


##########################################################################################################
#### Start the script
#########

### Set default command line options
retrieve_all=false
retrieve_counts=false
retrieve_tms=false
extract_flag=false
get_v8=false

### Parse command line options
usage="Retrieve and prepare GTEx RNA-seq data for in-house custom analyses.

USAGE: $(basename $0) [-h] [-a -c -t]
	where:
	-h Show this help text
	-a retrieve all expression data (mRNA counts and TPMs)
	-c retrieve only mRNA counts
	-t retrieve only mRNA TPMs
	-v get old V8 GTEx data used in the manuscript (default gets latest V10 data)

"
### parse all command line options
while getopts hactv opt; do
	case "$opt" in
		h) echo "$usage"; exit;;
		a) retrieve_all=true;;
		c) retrieve_counts=true;;
		t) retrieve_tms=true;;
		v) get_v8=true;;
		:)echo -e "Option -$OPTARG requires an argument.\n" >&2; echo "$usage" >&2; exit 1;;
		\?) echo ""; echo "$usage" >&2; exit 1;;
		*) echo ""; echo "$usage" >&2; exit 1;;
  esac
done

### if no expression set defined ask for one
if ! (${retrieve_all} || ${retrieve_counts} || ${retrieve_tms}); then 
	echo "No expression set defined. Use -a OR one of -c, -t" >&2
	echo""; echo "${usage}"
	exit 1
fi

cd ${DATA_FOLDER}
printf -v date '%(%Y-%m-%d)T' -1

gtex_web_base="https://storage.googleapis.com/adult-gtex"
if (${get_v8}); then
	GTEx_dir_id="GTEx-v8"
	### Add required variables for GTEx files V8
	sam_attr="GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
	web_sam_attr="${gtex_web_base}/annotations/v8/metadata-files/${sam_attr}"
	tpm_exp_data="GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct"
	web_tpm_exp_data="${gtex_web_base}/bulk-gex/v8/rna-seq/${tpm_exp_data}.gz"
	counts_exp_data="GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct"
	web_counts_exp_data="${gtex_web_base}/bulk-gex/v8/rna-seq/${counts_exp_data}.gz"
else
	GTEx_dir_id="GTEx-v10"
	### Add required variables for GTEx files V10
	sam_attr="GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"
	web_sam_attr="${gtex_web_base}/annotations/v10/metadata-files/${sam_attr}"
	tpm_exp_data="GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct"
	web_tpm_exp_data="${gtex_web_base}/bulk-gex/v10/rna-seq/${tpm_exp_data}.gz"
	counts_exp_data="GTEx_Analysis_v10_RNASeQCv2.4.2_gene_reads.gct"
	web_counts_exp_data="${gtex_web_base}/bulk-gex/v10/rna-seq/${counts_exp_data}.gz"
fi
mkdir -p ${GTEx_dir_id}
cd ${GTEx_dir_id}
working_folder=$(pwd)

### first get sample annotation file if not already present
echo "retrieving GTEx SampleAttributes data"
if test -f "${working_folder}/${sam_attr}"; then
	echo "${sam_attr} exists, moving on."
	echo ""
else
	wget ${web_sam_attr}
fi

#### iterate through the expression data types requested
if (${retrieve_all} || ${retrieve_tms}); then
	echo ""
	echo "retrieving GTEx TPM data"
	if test -f "${working_folder}/${tpm_exp_data}"; then
		echo "${tpm_exp_data} exists, moving on."
		echo ""
	else
		wget ${web_tpm_exp_data}
		extract_flag=true
	fi
	echo "Preparing for mRNA TPM data"
	process_folder ${date} ${sam_attr} ${tpm_exp_data} tpm
fi

if (${retrieve_all} || ${retrieve_counts}); then
	echo ""
	echo "retrieving GTEx count data"
	if test -f "${working_folder}/${counts_exp_data}"; then
		echo "${counts_exp_data} exists, moving on."
		echo ""
	else
		wget ${web_counts_exp_data}
		extract_flag=true
	fi
	echo "Preparing for mRNA Counts data"
	process_folder ${date} ${sam_attr} ${counts_exp_data} count
fi

### Sign off
echo ""
cd ${DATA_FOLDER}
echo "finished preparing all data for GTEx"
echo ""



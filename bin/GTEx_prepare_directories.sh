#!/usr/bin/env bash

##########################################################################################################
#### Prepare GTEx data for custom analyses
#### JToubia - January 2021
#########

##########################################################################################################
#### Source setup
#########
SCRIPT_FOLDER=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
BASE_FOLDER=`dirname ${SCRIPT_FOLDER}`
DATA_FOLDER="/data/projects/GTEx"
REF_FILES_FOLDER="${BASE_FOLDER}/REF_FILES"

### Add required variables for GTEx files
sam_attr="GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
tpm_exp_data="GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct"
counts_exp_data="GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct"

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

### Parse command line options
usage="Retrieve and prepare GTEx RNA-seq data for in-house custom analyses.

USAGE: $(basename $0) [-h] [-a -c -t]
	where:
	-h Show this help text
	-a retrieve all expression data (mRNA counts and TPMs)
	-c retrieve only mRNA counts
	-t retrieve only mRNA TPMs

"
### parse all command line options
while getopts hact opt; do
	case "$opt" in
		h) echo "$usage"; exit;;
		a) retrieve_all=true;;
		c) retrieve_counts=true;;
		t) retrieve_tms=true;;
		:)echo -e "Option -$OPTARG requires an argument.\n" >&2; echo "$usage" >&2; exit 1;;
		\?) echo ""; echo "$usage" >&2; exit 1;;
		*) echo ""; echo "$usage" >&2; exit 1;;
  esac
done

### if no expression set defined ask for one
if ! (${retrieve_all} || ${retrieve_counts} || ${retrieve_tms})
then 
	echo "No expression set defined. Use -a OR one of -c, -t" >&2
	echo""; echo "${usage}"
	exit 1
fi

cd ${DATA_FOLDER}
printf -v date '%(%Y-%m-%d)T' -1
GTEx_dir_id="GTEx-v8"
mkdir -p ${GTEx_dir_id}
cd ${GTEx_dir_id}
working_folder=`pwd`

### first get sample annotation file if not already present
echo "retrieving GTEx SampleAttributes data"
if test -f "${working_folder}/${sam_attr}"; then
	echo "${sam_attr} exists, moving on."
	echo ""
else
	wget https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt
fi

#### iterate through the expression data types requested
if (${retrieve_all} || ${retrieve_tms})
then
	echo ""
	echo "retrieving GTEx TPM data"
	if test -f "${working_folder}/${tpm_exp_data}"; then
		echo "${tpm_exp_data} exists, moving on."
		echo ""
	else
		wget https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz
		extract_flag=true
	fi
	echo "Preparing for mRNA TPM data"
	process_folder ${date} ${sam_attr} ${tpm_exp_data} tpm
fi

if (${retrieve_all} || ${retrieve_counts})
then
	echo ""
	echo "retrieving GTEx count data"
	if test -f "${working_folder}/${counts_exp_data}"; then
		echo "${counts_exp_data} exists, moving on."
		echo ""
	else
		wget https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz
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



#!/usr/bin/env bash

##########################################################################################################
#### Prepare RECOUNT3 data for custom analyses
#### JToubia - January 2023
#########

##########################################################################################################
#### Hardcoded DATA and SCRIPT dir loactions
#### Needs to be set on first run and should then remain in place
#### Modify "DATA_FOLDER" to suit your setup
#########
BASE_FOLDER="/home/jtoubia/Desktop/Projects/SRt"
DATA_FOLDER="${BASE_FOLDER}/RECOUNT3"
SCRIPT_FOLDER="${BASE_FOLDER}/bin/RECOUNT3_scripts"
REF_FILES_FOLDER="${BASE_FOLDER}/REF_FILES"

##########################################################################################################
#### Functions
#########
process_folder( ) {
	pID=$1
	date_stamp=$2
	
	create_ind_summary
}

create_ind_summary( ) {
	cd ${working_folder}
	Rscript --vanilla ${SCRIPT_FOLDER}/r3_prepare_quantification_file.R ${pID} ${date_stamp}
}


##########################################################################################################
#### Start the script
#########

### Set default command line options
pflag=false

### Parse command line options
usage="Retrieve and prepare RECOUNT3 RNA-seq data for in-house custom analyses.

USAGE: $(basename $0) [-h] -p <RECOUNT3 Project ID> 
	where:
	-h Show this help text
	-p RECOUNT3 project id: needs to be valid RECOUNT3 project id (ie. BRCA for TCGA data OR breast for GTEx). Required
	
To retrieve and prepare more than one RECOUNT dataset use a bash for loop like this;
for r3_code in COAD SARC LAML; do R3_prepare_directories.sh -a -p \$r3_code ; done

"
### parse all command line options
while getopts hp: opt; do
	case "$opt" in
		h) echo "$usage"; exit;;
		p) pflag=true; project_id=$OPTARG;;
		:)echo -e "Option -$OPTARG requires an argument.\n" >&2; echo "$usage" >&2; exit 1;;
		\?) echo ""; echo "$usage" >&2; exit 1;;
		*) echo ""; echo "$usage" >&2; exit 1;;
  esac
done

### check that -p has been used
if ! ${pflag}
then
	echo "-p must be used to pass a vaild RECOUNT3 Project ID" >&2
	echo""; echo "${usage}"
	exit 1
fi

### finally, check if user input a valid RECOUNT3 code
if grep -Fxq "${project_id}" <(cut -f 1 ${SCRIPT_FOLDER}/R3_API/R3_Study_Abbreviations.tsv)
then
	### Sign on
	echo "Starting retrieval of data and preparation of local repository for ${project_id}"
else
	echo "${project_id} not a valid RECOUNT3 study code, check your choice against the list below and try again"
	echo ""
	cat ${SCRIPT_FOLDER}/R3_API/R3_Study_Abbreviations.tsv
	exit 1
fi

cd ${DATA_FOLDER}
printf -v date '%(%Y-%m-%d)T' -1
r3_id_mod="${project_id}-${date}"

### setup parent directory for project
mkdir -p ${project_id}
cd ${project_id}
working_folder=`pwd`

### get the data
	echo ""
	process_folder ${project_id} ${date}

### warn the user about cached files
echo ""
echo "In the process of retrieving these files data were stored in cache for reuse."
echo "However, this can take up disc space which might be undesirable."
echo "To delete everything in this cache, start an R session and run the following command: recount3_cache_rm()"

### Sign off
echo ""
cd ${DATA_FOLDER}
echo "finished preparing all data for ${project_id}"
echo ""






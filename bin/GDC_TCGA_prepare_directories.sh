#!/usr/bin/env bash

##########################################################################################################
#### Prepare TCGA data for custom analyses
#### JToubia - January 2021
#########

##########################################################################################################
#### Source setup
#########
SCRIPT_FOLDER=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
BASE_FOLDER=`dirname ${SCRIPT_FOLDER}`
DATA_FOLDER="${BASE_FOLDER}/data/GDC"
REF_FILES_FOLDER="${BASE_FOLDER}/REF_FILES"

##########################################################################################################
#### Functions
#########
process_folder( ) {
	folder_name=$1
	gene_exp_folder_name=$2
	file_ext_orig=$3
	file_ext_new=$4
	tcga_id_mod_plus=$5
	rna_type=$6
	
#	echo "Processing ${folder_name}"
	get_files
	process_files
	#if [ ${rna_type} = "mRNA" ]; then extract_files; fi # no longer required
	rename_files
	get_clin_data
	create_summary
	if ! ${keep_data} ; then cleanup; fi
	cd ${working_folder}
}

get_files( ) {
	echo "getting files"
	mkdir -p ${folder_name}
	cd ${folder_name}
	${SCRIPT_FOLDER}/GDC_API/Download_SampleSheet.py ${project_id} ${folder_name}
	${SCRIPT_FOLDER}/GDC_API/Download_Files.py ${project_id} ${folder_name}
	tar -xif gdc_download_*.tar.gz
}

process_files( ) {
	echo "organising folders and processing files"
	mkdir -p .original
	mv * .original/
	mv .original original
	mkdir -p ${gene_exp_folder_name}
	cp original/*/*${file_ext_orig} ${gene_exp_folder_name}
	cd ${gene_exp_folder_name}
}

extract_files( ) {
	echo "unziping files"
	for file in ./*.gz; do 
		gunzip ${file}
	done
}

rename_files( ) {
	echo "renaming files"
	sample_tsv=(${working_folder}/${folder_name}/original/sample_sheet.tsv)
	${SCRIPT_FOLDER}/gdc_tcga_prepare_rename_files.py ${sample_tsv} ${file_ext_orig} ${file_ext_new}
}

get_clin_data( ) {
	echo "getting clinical data"
	cd ${working_folder}/${folder_name}
	mani_tsv=(${working_folder}/${folder_name}/original/MANIFEST.txt)
	${SCRIPT_FOLDER}/GDC_API/gdc-tsv-tool.py -c -o "${tcga_id_mod}_${rna_type}_clinical" ${mani_tsv}
}

create_summary( ) {
	echo "creating summary"
	cd ${working_folder}/${folder_name}
	${SCRIPT_FOLDER}/gdc_tcga_prepare_quantification_dir_to_summary.py \
		"${working_folder}/${folder_name}/${gene_exp_folder_name}/*${file_ext_new}" \
		${tcga_id_mod_plus} ${rna_type} ${REF_FILES_FOLDER} ${non_tcga_data}
}

combine_summary( ) {
	cd ${working_folder}
	${SCRIPT_FOLDER}/gdc_tcga_prepare_combine_summary.py \
	"${working_folder}/mRNA_expression_counts/GDC_${tcga_id_mod}_FPKM-mRNA_all.tsv" \
	"${working_folder}/isomiR_expression_rpm/GDC_${tcga_id_mod}_RPM-miRNAisoform_all.tsv" ${tcga_id_mod} 
}

cleanup( ) {
	echo "removing data that is no longer required"
	cd ${working_folder}/${folder_name}
	rm -rf original
	if [ -d "mRNA_gene_expr_files" ]
	then
		rm -rf mRNA_gene_expr_files
	else
		rm -rf miR_gene_expr_files
	fi
}


##########################################################################################################
#### Start the script
#########

### Set default command line options
pflag=false
retrieve_all=false
retrieve_counts=false
retrieve_miR_rpms=false
retrieve_isomiR_rpms=false
keep_data=false
non_tcga_data=false

### Parse command line options
usage="Retrieve and prepare TCGA RNA-seq data for in-house custom analyses.

USAGE: $(basename $0) [-h] -p <TCGA Project ID> [-a -c -r -R -k -n]
	where:
	-h Show this help text
	-p TCGA project id: needs to be valid TCGA project id as at the GDC (ie. TCGA-BRCA). Required
	-a retrieve all expression data (mRNA counts and TPMs as well as miR and isomiR RPMs)
	-c retrieve only mRNA counts
	-r retrieve only miR RPMs
	-R retrieve only isomiR RPMs
	-k keep all data (Default: False)
	-n data is not from TCGA study (Default: False)

To retrieve and prepare more than one TCGA cancer dataset use a bash for loop like this;
for tcga_code in TCGA-COAD TCGA-SARC TCGA-LAML; do GDC_TCGA_prepare_directories.sh -p \$tcga_code; done

To retrieve and prepare all TCGA cancer datasets you can loop through all lines in GDC_API/TCGA_Study_Abbreviations.tsv (WARNING: requires lots of time, network and disc space)
while read tcga_code; do GDC_TCGA_prepare_directories.sh -p \$tcga_code; done < <(cut -f1 GDC_API/TCGA_Study_Abbreviations.tsv | tail -n +2)
"
### parse all command line options
while getopts hp:acrRkn opt; do
	case "$opt" in
		h) echo "$usage"; exit;;
		p) pflag=true; project_id=$OPTARG;;
		a) retrieve_all=true;;
		c) retrieve_counts=true;;
		r) retrieve_miR_rpms=true;;
		R) retrieve_isomiR_rpms=true;;
		k) keep_data=true;;
		n) non_tcga_data=true;;
		:)echo -e "Option -$OPTARG requires an argument.\n" >&2; echo "$usage" >&2; exit 1;;
		\?) echo ""; echo "$usage" >&2; exit 1;;
		*) echo ""; echo "$usage" >&2; exit 1;;
  esac
done

### check that -p has been used
if ! ${pflag}
then
	echo "-p must be used to pass a vaild TCGA Project ID" >&2
	echo""; echo "${usage}"
	exit 1
fi

### if no expression set defined ask for one
if ! (${retrieve_all} || ${retrieve_counts} || ${retrieve_miR_rpms} || ${retrieve_isomiR_rpms})
then 
	echo "No expression set defined. Use -a OR one or more of -c, -r, -R" >&2
	echo""; echo "${usage}"
	exit 1
fi

### finally, check if user input a valid TCGA/GDC code
if grep -Fxq "${project_id}" <(cut -f 1 ${REF_FILES_FOLDER}/study_abbreviations/TCGA_Study_Abbreviations.tsv)
then
	### Sign on
	echo "Starting retrieval of data and preparation of local repository for ${project_id}"
else
	echo "${project_id} not a valid TCGA study code, check your choice against the list below and try again"
	echo ""
	cat ${REF_FILES_FOLDER}/study_abbreviations/TCGA_Study_Abbreviations.tsv
	exit 1
fi

cd ${DATA_FOLDER}
printf -v date '%(%Y-%m-%d)T' -1
tcga_id_mod="${project_id}-${date}"

### setup parent directory for project
mkdir -p ${project_id#*TCGA-}
cd ${project_id#*TCGA-}
working_folder=`pwd`

#### iterate through the expression data types requested
if (${retrieve_all} || ${retrieve_counts})
then
	echo ""
	echo "Preparing for mRNA Counts data"
	process_folder mRNA_expression_counts mRNA_gene_expr_files _star_gene_counts.tsv .star_gene_counts.tsv ${tcga_id_mod}_Count-mRNA mRNA
fi

if (${retrieve_all} || ${retrieve_miR_rpms})
then
	echo ""
	echo "Preparing for miR RPM data"
	process_folder miR_expression_rpm miR_gene_expr_files .quantification.txt .quantification.txt ${tcga_id_mod}_RPM-miRNA miRNA
fi

if (${retrieve_all} || ${retrieve_isomiR_rpms})
then
	echo ""
	echo "Preparing for isomiR RPM data"
	process_folder isomiR_expression_rpm miR_gene_expr_files .quantification.txt .quantification.txt ${tcga_id_mod}_RPM-miRNAisoform isomiRNA
fi

### Finally, create a combined mRNA-miRNA spreadsheet
if (${retrieve_all}) || (${retrieve_counts} && ${retrieve_isomiR_rpms})
then
	echo ""
	echo "Combining mRNA FPKM and isomiR RPM data"
	combine_summary
fi

### Sign off
echo ""
cd ${DATA_FOLDER}
echo "finished preparing all data for ${project_id}"
echo ""






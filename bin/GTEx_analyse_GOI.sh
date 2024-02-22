#!/usr/bin/env bash

##########################################################################################################
#### DE analysis by GOI stratification in GTEx data
#### JToubia - January 2021
#########

##########################################################################################################
#### Source setup
#########
SCRIPT_FOLDER=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
BASE_FOLDER=`dirname ${SCRIPT_FOLDER}`
DATA_FOLDER="/data/projects/GTEx"
REF_FILES_FOLDER="${BASE_FOLDER}/REF_FILES"
GSEA_EXE="${SCRIPT_FOLDER}/GSEA_Linux_4.2.3/gsea-cli.sh"

### and record where we start from
#working_folder=`pwd`

##########################################################################################################
#### Functions
#########
corr_analysis( ) {
	echo "Preparing to run Corr Analysis for ${GOI}"
	GTEx_z2n=(${DATA_FOLDER}/${GTEx_code}/GTEx-*_tpm-mRNA.tsv)
	Rscript --vanilla ${SCRIPT_FOLDER}/gtex_analyse_correlation_plots.R ${GOI} ${output_folder} ${GTEx_z2n[0]} ${REF_FILES_FOLDER}
	return_code=$?
	if [[ $return_code -ne 0 ]]
	then
		return $return_code
	fi 
	cd ${output_folder}
}

de_analysis( ) {
	echo "Preparing to run Stratification followed by DE Analysis for ${GOI}"
	GTEx_tpm=(${DATA_FOLDER}/${GTEx_code}/GTEx-*_tpm-mRNA.tsv)
	GTEx_counts=(${DATA_FOLDER}/${GTEx_code}/GTEx-*_count-mRNA.tsv)
	GTEx_rpm=""
	Rscript --vanilla ${SCRIPT_FOLDER}/gtex_analyse_stratify_edgeR.R ${GOI} ${output_folder} ${strat_by} ${percentile} ${REF_FILES_FOLDER} ${GSEA_EXE} ${GTEx_tpm[0]} ${GTEx_counts[0]} ${GTEx_rpm[0]}
	return_code=$?
	if [[ $return_code -ne 0 ]]
	then
		return $return_code
	fi 
	cd ${output_folder}
}


##########################################################################################################
#### Start the process
#########

### Set default command line options
pflag=false
gflag=false
sflag=false
tflag=false
do_all=false
do_corr=false
do_de=false
output_flag=false

### Parse command line options
usage="Differential expression analysis of GTEx data stratified into high and low groups by gene of interest 
Please run this wrapper script in the directory of the desired output location

NOTE: Composite analyses can be run using the plus charater (+) for additive combinations or
by using the modulus character (%) for ratio. The two special characters are used between gene names like so
Additive example: ESR1+PGR+ERBB2 or Ratio example: ESRP1%ZEB1

USAGE: $(basename $0) [-h] -p <GTExTissueID> -g <GOI> -s <StratifyBy> -t <Percentile> -a -c -d -o
	where:
	-h Show this help text
	-p GTEx tissue id: needs to be valid GTEx tissue id as at GTEx (ie. Breast). Required
	-g Gene of interest: needs to be a valid HGNC symbol (ie. ZEB1). Required
	-s Stratify by molecule: Must match -g and can only be mRNA at present. Required
	-t Percentile: startify data into top and bottom x percentil (valid x between 2 and 25). Required
	-a Do all analyses
	-c Do correlation analysis only
	-d Do differential expression analysis only 
	-o Path to where the analysis should be stored
"

### parse all command line options
while getopts hg:p:s:t:o:acd opt; do
	case "$opt" in
		h) echo "$usage"; exit;;
		p) pflag=true; project_id=$OPTARG;;
		g) gflag=true; GOI=$OPTARG;;
		s) sflag=true; strat_by=$OPTARG;;
		t) tflag=true; percentile=$OPTARG;;
		o) output_flag=true; output_folder=$OPTARG;;
		a) do_all=true;;
		c) do_corr=true;;
		d) do_de=true;;
		:)echo -e "Option -$OPTARG requires an argument.\n" >&2; echo "$usage" >&2; exit 1;;
		\?) echo ""; echo "$usage" >&2; exit 1;;
		*) echo ""; echo "$usage" >&2; exit 1;;
  esac
done

### check that all required flags have been given
#if ! (${pflag} && ${gflag} && ${sflag} && ${tflag} && ${output_flag})
#then
#	echo "Missing arguments. Please check your commanline call" >&2
#	echo""; echo "${usage}"
#	exit 1
#fi

### check that all required flags have been given
# changed such that we no longer need to provide s/tflag for correlation analysis only
if (${do_all} || ${do_de})
then 
	if ! (${pflag} && ${gflag} && ${sflag} && ${tflag} && ${output_flag})
	then
		echo "Missing arguments. Please check your commanline call" >&2
		echo""; echo "${usage}"
		exit 1
	fi
else
	if ! (${pflag} && ${gflag} && ${output_flag})
	then
		echo "Missing arguments. Please check your commanline call" >&2
		echo""; echo "${usage}"
		exit 1
	fi
fi

### check if user input a valid GTEx code
if grep -Fxq "${project_id}" <(cut -f 1 ${REF_FILES_FOLDER}/study_abbreviations/GTEx_Study_Abbreviations.tsv)
then
	GTEx_code=${project_id#*GTEx-}
else
	echo "${project_id} not a valid GTEx study code, check your choice against the list below and try again"
	echo ""
	cat ${REF_FILES_FOLDER}/study_abbreviations/GTEx_Study_Abbreviations.tsv
	exit 1
fi

### check if user input a valid gene name
if [[ "${GOI}" == *"+"* ]] || [[ "${GOI}" == *"%"* ]]
then
	echo "complex analysis selected, not checking gene names yet"
else
	if grep -Fxq "${GOI^^}" <(cut -f 2 ${REF_FILES_FOLDER}/gencode.v26.annotation.lookup)
	then
		GOI=${GOI^^}
	else
		echo "${GOI^^} is not a valid gene name or is not in the current build, check your spelling maybe?"
		echo "ALternatively, look in this file (${REF_FILES_FOLDER}/gencode.v26.annotation.lookup)"
		echo "for a complete list of valid gene names"
		exit 1
	fi
fi

### check if user input a valid RNA species
if (${do_all} || ${do_de})
then
	if [ "${strat_by}" != "mRNA" ]
	then
		echo "Stratification can only be done by mRNA" >&2
		echo""; echo "${usage}"
		exit 1
	fi

	### check if user input a valid integer
	if ((${percentile} < 2 || ${percentile} > 25))
	then
		echo "Percentile must be an integer in the range 2 - 25" >&2
		echo""; echo "${usage}"
		exit 1
	fi
fi

### if no analysis defined ask for one
if ! (${do_all} || ${do_corr} || ${do_de})
then 
	echo "No analysis defined. Use -a OR one or more of -c, -d" >&2
	echo""; echo "${usage}"
	exit 1
fi

### check if user input a complex gene name (additive or ratio)
if [[ "${GOI}" == *"+"* ]] || [[ "${GOI}" == *"%"* ]]
then
	echo "Looks like a complex expression analysis, deactivating Corr Analyses" >&2
	do_all=false
	do_corr=false
	do_de=true
fi

### Sign on
echo "Starting Correlation and DE analysis for ${GOI} using GTEx subset: ${project_id}"

### iterate through the objectives
if (${do_all} || ${do_corr})
then
	echo "Corr Analysis selected"
	corr_analysis
	if [[ $return_code -ne 0 ]]
	then
		exit $return_code
	fi 
fi

if (${do_all} || ${do_de})
then
	echo "DE Analysis selected"
	de_analysis
	if [[ $return_code -ne 0 ]]
	then
		exit $return_code
	fi 
fi

### Final clean up
cd $output_folder
find . -type d -empty -delete

### Sign off
echo "Finished running all analyses for ${GOI}"









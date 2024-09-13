#!/usr/bin/env bash

##########################################################################################################
#### DE analysis by GOI stratification in RECOUNT3 data
#### JToubia - January 2023
#########

##########################################################################################################
#### Source setup
#########
SCRIPT_FOLDER=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
BASE_FOLDER=`dirname ${SCRIPT_FOLDER}`
DATA_FOLDER="${BASE_FOLDER}/data/RECOUNT3"
REF_FILES_FOLDER="${BASE_FOLDER}/REF_FILES"
GSEA_EXE="${SCRIPT_FOLDER}/GSEA/gsea-cli.sh"

### and record where we start from
working_folder=`pwd`

##########################################################################################################
#### Functions
#########
corr_analysis( ) {
	echo "Preparing to run Corr Analysis for ${GOI}"
	R3_z2n=(${DATA_FOLDER}/${R3_code}/R3-*_tpm-mRNA.tsv)
	Rscript --vanilla ${SCRIPT_FOLDER}/r3_analyse_correlation_plots.R ${GOI} ${R3_z2n[0]} ${REF_FILES_FOLDER}
	cd ${working_folder}
}

de_analysis( ) {
	echo "Preparing to run Stratification followed by DE Analysis for ${GOI}"
	R3_tpm=(${DATA_FOLDER}/${R3_code}/R3-*_tpm-mRNA.tsv)
	R3_counts=(${DATA_FOLDER}/${R3_code}/R3-*_count-mRNA.tsv)
	R3_rpm=""
	Rscript --vanilla ${SCRIPT_FOLDER}/r3_analyse_stratify_edgeR.R ${GOI} ${strat_by} ${percentile} ${switch} ${REF_FILES_FOLDER} ${GSEA_EXE} ${R3_tpm[0]} ${R3_counts[0]} ${R3_rpm[0]}
	cd ${working_folder}
}

trace () {
stamp=`date +%Y-%m-%d_%H:%M:%S`
echo $stamp: $* >> SCA_command.log
}

##########################################################################################################
#### Start the process
#########
trace Starting SCA
trace $0 $@

### Set default command line options
pflag=false
gflag=false
sflag=false
tflag=false
gsea=false
switch=false
do_all=false
do_corr=false
do_de=false

### Parse command line options
usage="Differential expression analysis of RECOUNT3 data stratified into high and low groups by gene of interest 
Please run this wrapper script in the directory of the desired output location

NOTE: Composite analyses can be run using the plus charater (+) for additive combinations or
by using the modulus character (%) for ratio. The two special characters are used between gene names like so
Additive example: ESR1+PGR+ERBB2 or Ratio example: ESRP1%ZEB1

USAGE: $(basename $0) [-h] -p <RECOUNT3 projectID> -g <GOI> -s <StratifyBy> -t <Percentile> -e -S -a -c -d
	where:
	-h Show this help text
	-p RECOUNT3 tissue id: needs to be valid RECOUNT3 tissue id as at RECOUNT3 (ie. BRCA for TCGA or BREAST for GTEx). Required
	-g Gene of interest: needs to be a valid HGNC symbol (ie. ZEB1). Required
	-s Stratify by molecule: Must match -g and can only be mRNA at present. Required
	-t Percentile: startify data into top and bottom x percentile (valid x between 2 and 25). Required
	-e Enrichment analyses: Run GSEA on DE results (Default: Only run WebGestalt)
	-S Switch pairwise comparison: find genes DE in low group compared to high group (Default: high compared to low)
	-a Do all analyses
	-c Do correlation analysis only
	-d Do differential expression analysis only 
"

### parse all command line options
while getopts hg:p:s:t:eSacd opt; do
	case "$opt" in
		h) echo "$usage"; exit;;
		p) pflag=true; project_id=$OPTARG;;
		g) gflag=true; GOI=$OPTARG;;
		s) sflag=true; strat_by=$OPTARG;;
		t) tflag=true; percentile=$OPTARG;;
		e) gsea=true;;
		S) switch=true;;
		a) do_all=true;;
		c) do_corr=true;;
		d) do_de=true;;
		:)echo -e "Option -$OPTARG requires an argument.\n" >&2; echo "$usage" >&2; exit 1;;
		\?) echo ""; echo "$usage" >&2; exit 1;;
		*) echo ""; echo "$usage" >&2; exit 1;;
  esac
done

### check that all required flags have been given
if ! (${pflag} && ${gflag} && ${sflag} && ${tflag})
then
	echo "Missing arguments. Please check your commanline call" >&2
	echo""; echo "${usage}"
	exit 1
fi

### check if user input a valid RECOUNT3 code
if grep -Fxq "${project_id}" <(cut -f 1 ${REF_FILES_FOLDER}/study_abbreviations/R3_Study_Abbreviations.tsv)
then
	R3_code=${project_id}
else
	echo "${project_id} not a valid RECOUNT3 study code, check your choice against the list below and try again"
	echo ""
	cat ${REF_FILES_FOLDER}/study_abbreviations/R3_Study_Abbreviations.tsv
	exit 1
fi

### check if user input a valid gene name
if [[ "${GOI}" == *"+"* ]] || [[ "${GOI}" == *"%"* ]]
then
	echo "Complex analysis selected, not checking gene names, deactivating Corr Analyses"
	do_all=false
	do_corr=false
	do_de=true
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
if [ "${strat_by}" != "mRNA" ]
then
	echo "Stratification can only be done by mRNA" >&2
	echo""; echo "${usage}"
	exit 1
fi

### check if user input a valid integer
if ((${percentile%.*} < 2 || ${percentile%.*} > 25))
then
	echo "Percentile must be an integer in the range 2 - 25" >&2
	echo""; echo "${usage}"
	exit 1
fi

### if no analysis defined ask for one
if ! (${do_all} || ${do_corr} || ${do_de})
then 
	echo "No analysis defined. Use -a OR one or more of -c, -d" >&2
	echo""; echo "${usage}"
	exit 1
fi

### check if enrichment analysis selected otherwise clear GSEA_EXE variable
if (${gsea})
then 
	echo "Enrichment analyses requested. Note: run time will be extended significantly"
else
	GSEA_EXE="false"
fi

### Sign on
echo "Starting analysis for ${GOI} using RECOUNT3 subset: ${project_id}"

### iterate through the objectives
if (${do_all} || ${do_corr})
then
	echo "Corr Analysis selected"
	corr_analysis
fi

if (${do_all} || ${do_de})
then
	echo "DE Analysis selected"
	
	### check if user requested to switch DE groups
	if (${switch})
	then 
		echo "Switch groups requested, evaluating changes in low group compared to high"
	fi
	
	de_analysis
fi

### Final clean up
cd $working_folder
find . -type d -empty -delete
if (${do_all} || ${do_de})
then
	cd DE_Analysis
	${SCRIPT_FOLDER}/post_analysis_organisation.sh
	cd $working_folder
fi

trace Finished SCA

### Sign off
echo "Finished running all analyses for ${GOI}"









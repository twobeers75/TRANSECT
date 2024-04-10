#!/usr/bin/env bash

##########################################################################################################
#### Script to test Prepare and Analyse parts of the pipeline with different parameters  
#### as well as reproducing the results in the manuscript
#########

#### MODIFY THE FOLLOWING PATH TO SUIT YOUR SETUP
SCA_HOME="/data_b/tools/SCA"

#### Should not require any modifications below this line
SCA_bindir="${SCA_HOME}/bin"
SCA_outdir="${SCA_HOME}/output"

##########################################################################################################
#### Get data - If you already have these you can comment out the following 3 commands
#########
# RECOUNT3
for r3_code in BRCA BREAST PRAD PROSTATE LAML BLOOD; do ${SCA_bindir}/R3_prepare_directories.sh -p $r3_code; done
# GDC-TCGA
for tcga_code in TCGA-BRCA TCGA-PRAD TCGA-LAML; do ${SCA_bindir}/GDC_TCGA_prepare_directories.sh -p $tcga_code -a; done
# GTEx
${SCA_bindir}/GTEx_prepare_directories.sh -a; done

##########################################################################################################
#### Single test
#########
test_prostate="${SCA_outdir}/TestManData_OG2/Prostate"
mkdir -p ${test_prostate}
cd ${test_prostate}
mkdir -p R3_PRAD_ZEB1 R3_PROSTATE_ZEB1 TCGA_PRAD_ZEB1 GTEx_Prostrate_ZEB1

# RECOUNT3 Single ZEB1 with corr analysis
cd ${test_prostate}/R3_PRAD_ZEB1
${SCA_bindir}/R3_analyse_GOI.sh -p PRAD -g ZEB1 -s mRNA -t 5 -a -e
echo ""; echo ""
cd ${test_prostate}/R3_PROSTATE_ZEB1
${SCA_bindir}/R3_analyse_GOI.sh -p PROSTATE -g ZEB1 -s mRNA -t 10 -a -e
echo ""; echo ""
# GDC Single ZEB1 with corr analysis
cd ${test_prostate}/TCGA_PRAD_ZEB1
${SCA_bindir}/GDC_TCGA_analyse_GOI.sh -p TCGA-PRAD -g ZEB1 -s mRNA -t 5 -a -e
echo ""; echo ""
# GTEx Single ZEB1 with corr analysis
cd ${test_prostate}/GTEx_Prostrate_ZEB1
${SCA_bindir}/GTEx_analyse_GOI.sh -p Prostate -g ZEB1 -s mRNA -t 10 -a -e
echo ""; echo ""

##########################################################################################################
#### Additive test
#########
test_breast="${SCA_outdir}/TestManData_OG2/Breast"
mkdir -p ${test_breast}
cd ${test_breast}
mkdir -p R3_BRCA_ESR1+PGR+ERBB2 R3_BREAST_ESR1+PGR+ERBB2 TCGA_BRCA_ESR1+PGR+ERBB2 GTEx_Breast_ESR1+PGR+ERBB2

# RECOUNT3 Additive TNBC genes with switch
cd ${test_breast}/R3_BRCA_ESR1+PGR+ERBB2
${SCA_bindir}/R3_analyse_GOI.sh -p BRCA -g ESR1+PGR+ERBB2 -s mRNA -t 3 -S -d -e
echo ""; echo ""
cd ${test_breast}/R3_BREAST_ESR1+PGR+ERBB2
${SCA_bindir}/R3_analyse_GOI.sh -p BREAST -g ESR1+PGR+ERBB2 -s mRNA -t 6 -S -d -e
echo ""; echo ""
# GDC Additive TNBC genes with switch
cd ${test_breast}/TCGA_BRCA_ESR1+PGR+ERBB2
${SCA_bindir}/GDC_TCGA_analyse_GOI.sh -p TCGA-BRCA -g ESR1+PGR+ERBB2 -s mRNA -t 3 -S -d -e
echo ""; echo ""
# GTEx Additive TNBC genes with switch
cd ${test_breast}/GTEx_Breast_ESR1+PGR+ERBB2
${SCA_bindir}/GTEx_analyse_GOI.sh -p Breast -g ESR1+PGR+ERBB2 -s mRNA -t 6 -S -d -e
echo ""; echo ""

# additionally do the ER+PGR- analysis
cd ${test_breast}
mkdir -p R3_BRCA_ESR1-PGR
cd R3_BRCA_ESR1-PGR
${SCA_bindir}/R3_analyse_GOI.sh -p BRCA -g ESR1%PGR -s mRNA -t 2 -d -e
echo ""; echo ""

##########################################################################################################
#### Ratio test
#########
test_blood="${SCA_outdir}/TestManData_OG2/Blood"
mkdir -p ${test_blood}
cd ${test_blood}
mkdir -p R3_LAML_IL3RA-CSF2RB R3_BLOOD_IL3RA-CSF2RB TCGA_LAML_IL3RA-CSF2RB GTEx_Blood_IL3RA-CSF2RB

# RECOUNT3 Additive TNBC genes with switch
cd ${test_blood}/R3_LAML_IL3RA-CSF2RB
${SCA_bindir}/R3_analyse_GOI.sh -p LAML -g IL3RA%CSF2RB -s mRNA -t 12 -d -e
echo ""; echo ""
cd ${test_blood}/R3_BLOOD_IL3RA-CSF2RB
${SCA_bindir}/R3_analyse_GOI.sh -p BLOOD -g IL3RA%CSF2RB -s mRNA -t 3 -d -e
echo ""; echo ""
# GDC Additive TNBC genes with switch
cd ${test_blood}/TCGA_LAML_IL3RA-CSF2RB
${SCA_bindir}/GDC_TCGA_analyse_GOI.sh -p TCGA-LAML -g IL3RA%CSF2RB -s mRNA -t 12 -d -e
echo ""; echo ""
# GTEx Additive TNBC genes with switch
cd ${test_blood}/GTEx_Blood_IL3RA-CSF2RB
${SCA_bindir}/GTEx_analyse_GOI.sh -p Blood -g IL3RA%CSF2RB -s mRNA -t 3 -d -e
echo ""; echo ""

















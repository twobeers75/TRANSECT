<img title="TRANSECT" alt="TRANSECT" src="https://github.com/twobeers75/TRANSECT/blob/main/TRANSECT_t.png">

#
| [Installation](#installation) | [Usage](#usage) | [Output](#output) | [Manual](#manual) |


## TRANSECT Introduction

TRANSECT stratifies cohort data into user defined strata (groups) based solely on expression, and subsequently compares the strata one to the other. This version of TRANSECT is tuned specifically for the interrogation of gene/s from publicly available tissue-specific gene expression data.

Currently, there are many global efforts to collect and collate tissue-specific gene expression (RNA-sequencing) data from large cohorts of diseased and non-diseased subjects for public use (GTEx, TCGA, …, to name only a few). These data sets are used to bolster and aid in the investigation of gene specific expression and regulation as well as several other biological queries. One such application that is currently underutilised is the stratification and subsequent differential expression of single or composite genes in order to assess "transcriptome state". This type of analysis has only in recent times become feasible thanks to the large number of participants in cohort studies and, will only intensify in power and usefulness as these collections continue to grow in number and diversity. The results from this type of analysis can be used to rapidly investigate transcriptome state within natural physiological expression levels without economic burden. This can be applied before laboratory investigations and/or between the transition from *in vitro* to *in vivo* studies in order to assess the feasibility of further, possibly expensive, experimentation.

## TRANSECT Web-application
A web-accessible version is freely available for use [here](https://transect.au)

## Version

24.03

## Authors

John Toubia

## Hardware requirements

64 bit Linux (developed solely on Ubuntu 22.04, tested on 22.04 and 24.04)

10GB RAM (12+ recommended)

## Software requirements

- Python 3
  - see pip_requirements.txt for a complete list of required python packages
- R
  - see r_requirements.txt for a complete list of required R packages
- Java
  - for Gene Set Enrichment Analysis


## Installation
<span id="#installation"></span>

NOTE: TRANSECT requires and depends on numerous packages and applications. For this reason it is recommended for most use cases, to install TRANSECT within a Conda environment. This is not only easier and cleaner but also saves a significant amount of time (*from approx. 60 minutes to less than 10*). For those not wanting to use Conda, please see the instructions for an alternate installation method in the TRANSECT manual.

To start, clone the repo

```sh
git clone https://github.com/twobeers75/TRANSECT.git
```

Or, download the package and extract it

Browse to https://github.com/twobeers75/TRANSECT, click on the green "Code" button followed by "Download ZIP" (note the download location). 
Find the downloaded ZIP file and move it to an appropriate location if required before extracting the contents and renaming the folder

```sh
unzip TRANSECT-main.zip
mv TRANSECT-main TRANSECT
```

Install Conda (requires Conda version > 24.1.0). You can skip this step if you already have it installed on your system.

There are many wikis on how to install Conda for Ubuntu, [here](https://docs.anaconda.com/miniconda/) is just one.

*(approx. 1-2min)*

```sh
### follow the instructions outlined in the link above which should look something like this for Linux
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh

### don't forget to initialize the bash shell. Mac users need to check their default shell and change the following command appropriately
~/miniconda3/bin/conda init bash

*** Afterwards, you will be asked to restart your terminal whereby you should see (base) at the prompt. Ignore for now.
```

Create the TRANSECT Conda environment. In this step, we will create an environment with all the tools and dependencies required to run TRANSECT. 

*(approx. 5-10mins)*

```sh
### change into the top directory of the downloaded folder (TRANSECT) and navigate to the INSTALL folder
cd <path to>/TRANSECT/INSTALL

### First, run the conda install script (you may need to run "chmod 755 TRANSECT_conda_install.sh")
./TRANSECT_conda_install.sh

### Next, upon succesful completion of the previous step, activate the newly created environment
conda activate TRANSECT

### Finally whilst still within the TRANSECT/INSTALL directory, in the TRANSECT environment run the post installation script to complete the setup
./TRANSECT_post_conda_install.sh

# Note: you need to reactivate the TRANSECT environment at this point.
conda deactivate
conda activate TRANSECT
```

**And that's it!** You should now have all the necessary applications and dependencies in the TRANSECT environment to run this application. Please note, just like any virtual environment you are required to activate the TRANSECT environment in order to use the application. You can deactivate at will when not in use. 

A few extra commands for those not accustomed to Conda environments

```sh
### By default, Conda auto-activates the "base" default environment. 
### Each time you open a terminal you will automatically be within the (base) environment.
### I prefer not to have this happen. To disable, run the following commands
# first, deactivate any environment until there is no () at the begining of the prompt, then turn off auto_activate_base
conda deactivate
conda config --set auto_activate_base false

### You need to activate the TRANSECT environment each time you run TRANSECT
# to check what environments you have on your system
conda env list
# You should see "base" and "TRANSECT" at the very least

### To activate TRANSECT
conda activate TRANSECT

### Once you have finished with TRANSECT, deactivate the environment
conda deactivate
```

More info about managing Conda environment can be found [here](https://docs.conda.io/projects/conda/en/4.6.0/user-guide/tasks/manage-environments.html)



## Usage

<span id="#usage"></span> 

TRANSECT has two main operations; **Prepare** and **Analyse**. 

**Prepare** is a process that retrieves the raw data from online repositories and prepares it (if required) for analysis. TRANSECT comes bundled with three different prepare scripts, one each for RECOUNT3, GDC-TCGA and GTEx data. 

Example prepare command for RECOUNT3 TCGA PRAD cohort;

```sh
### First, if not already make sure to activate the TRANSECT environment
conda activate TRANSECT

### RECOUNT3 (approx. 5mins)
# change into the top directory of TRANSECT
cd <path to>/TRANSECT
# run the RECOUNT3 prepare script for TCGA-PRAD
bin/R3_prepare_directories.sh -p PRAD

### NOTE: you can use the -h parameter to see the full help menu. ie.
bin/R3_prepare_directories.sh -h

### The prepare commands for GDC-TCGA and GTEx are similar but not identical. Use -h and see the manual for full details.
```

Be aware that some of these collections are large and require substantial disk space. They will take a considerable amount of time to download and process too. For example, downloading and processing GDC TCGA-BRCA takes just over 30 minutes (using a high speed network connection and an up to date workstation) and requires more than 14GB of disk space (most of which can and by default is, deleted afterwards). In comparison, GDC TCGA-LAML takes less than 5 minutes to retrieve and less than 2GB of disc space. 

Preparation for both GDC-TCGA and the RECOUNT3 data is done individually by tissue type but can also be done in batch mode (see the relevant script help menu for instructions). Preparing GTEx data on the other hand, retrieves in bulk all tissue types in a single table before separating them into individual files based on tissue type (again, see the help menu for more details). 

All downloaded data is stored in "TRANSECT/data/<GDC|GTEx|RECOUNT3>" in appropriate folders.



**Analyse** is a process that uses the prepared public data from above conducts the stratified differential expression and produces all the outputs. Like with the prepare operation, TRANSECT comes bundled with three analyse scripts, one each for RECOUNT3, GDC-TCGA and GTEx data.

Example analyse command for the gene ZEB1 in the RECOUNT3 TCGA PRAD cohort;

```sh
### You should still be in the TRANSECT environment but if not, make sure to activate it again

### RECOUNT3 (approx. 5-10mins)
# change into the top directory of TRANSECT if not already there
# TRANSECT saves output in the current working folder so best to create a folder specifically for each run 
cd <path to>/TRANSECT/output/RECOUNT3
mkdir -p ZEB1_PRAD_test
cd ZEB1_PRAD_test
# Now, run the RECOUNT3 prepare script using the PRAD data we just retreived, investigating the gene ZEB1, with all outputs.
# because we are not in the script folder and have not added these to our PATH variable, you will need to provide the full path to the script
<full path to>/TRANSECT/bin/R3_analyse_GOI.sh -p PRAD -g ZEB1 -s mRNA -t 5 -a

### NOTE: use the -h parameter to see the full help menu. ie.
<full path to>/TRANSECT/bin/R3_analyse_GOI.sh -h

### The analysis commands for GDC-TCGA and GTEx are similar but not identical. Again, use -h or see the manual for full details.
```

## Output

<span id="#output"></span>

TRANSECT takes in a cohort dataset and processes the data as follows. 
1)	First, TRANSECT partition the data by the expression of a gene/s of interest into low and high strata 
2)	Subsequently, TRANSECT compares the resulting strata, one to the other, to identify differentially expressed genes
3)	And finally, TRANSECT uses the results from the DE analysis to run functional annotation and enrichment analyses

The outputs from TRANSECT are likewise grouped into 3 categories and returned in three folders in the working directory from where the program is executed

**01-Stratification**
The stratification process produces 2 tables, and 3 plots. 

1. GOI_exp_raw_OG.tsv contain the raw original expression (TPM) data for all gene/s of interest 
2. GOI_exp_with_strat.tsv contains the same data sorted with additional columns relating to the participants ranking score, percentiles and quantile values.
3. TPM_histogram.html or TPM_Boxplot_Sina.html or TPM_Scatter.html, all which plot data from the two tables above differently depending on the chosen TRANSECT mode in an attempt to describe the distribution of gene expression across the cohort participants
4. TPM_N-T_boxplot.html which shows the distribution of expression partitioned by disease state when available 
5. TPM_strat_boxplot.html which plots the low and high strata participants resulting from the stratification process

**02-DE**
The DE analysis produces many tables and plots most easily described as follows.

1. DE Setup – design.tsv and gene_raw_expression_data_cpm.csv 
2. DE QC – bcv and mean_var.png plots as well as the MDS-Plot.html in the glimma-plots folder
3. Normalised expression tables - gene_normalised_expression_data_cpm.csv (also in log form)
4. DE result tables - High_Vs_Low _de_sigFC.csv and top_tags.csv
5. DE result plots - High_Vs_Low_volcano.png and High_Vs_Low_heatmap.png as well as an interactive version of the volcano plot in the glimma-plots folder
6. The glimma-plots folder containing the interactive web plots and associated data

**03-Enrichment**
The 2 enrichment analyses result in the production of two folders each with a separate collection of tables and plots.

1. GSEA

   When selected, this folder contains the output folders from running GSEA against the Hallmark as well as the Curated MSigDB collections respectively. Within each folder, users can open the index.html file to access and interact with the results in a web browser. In addition, the results are summarised and provided in tabular form (.csv) as well as interactive form (.html)

   GSEA input data – 3 text files used for the GSEA analysis are saved in the top-level folder. The default GSEA method used by TRANSECT is the pre-ranked method. Input for this analysis can be found in the .rnk file. Provided but not used by TRANSECT are alternate GSEA input files (.cls and .txt).

2. WebGestalt

   The ORA results are presented in six folders; two each for disease, gene ontology and pathway enrichment, for up and down regulated genes separately (when available). Within each folder, users can open the .html file to access and interact with the results in a web browser.

## Manual

<span id="#manual"></span>

...

## Publication

Citation details to come (hopefully!)

## Licence

TRANSECT is made available under the terms of the MIT license, a copy of which is included in the distribution in the LICENSE file.

GSEA is made available under the terms of a BSD-style license, a copy of which is included in the GSEA folder in the LICENSE.txt file. See that file for exact terms and conditions.

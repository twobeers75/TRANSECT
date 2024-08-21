<img title="TRANSECT" alt="TRANSECT" src="https://github.com/twobeers75/TRANSECT/blob/main/TRANSECT_t.png">

#
| [Installation](#installation) | [Usage](#usage) | [Output](#output) | [Manual](#manual) |


## TRANSECT Introduction

TRANSECT is a bioinformatics application designed to run in a Linux Terminal. In brief, TRANSECT stratifies large cohort data into user defined subsets based solely on expression, which are subsequently compared one to the other. This version of TRANSECT is tuned specifically for the interrogation of gene/s from publicly available tissue-specific gene expression data.

Currently, there are many global efforts to collect and collate tissue-specific gene expression (RNA-sequencing) data from large cohorts of diseased and non-diseased subjects for public use (GTEx, TCGA, …, to name only a few). These data sets are used to bolster and aid in the investigation of gene specific expression and regulation as well as several other biological queries. One such application that is currently underutilised is the stratification and subsequent differential expression of single or composite genes in order to assess "transcriptome state". This type of analysis has only in recent times become feasible thanks to the large number of participants in cohort studies and, will only intensify in power and usefulness as these collections continue to grow in number and diversity. The results from this type of analysis can be used to rapidly investigate transcriptome state within natural physiological expression levels without economic burden. This can be applied before laboratory investigations and/or between the transition from *in vitro* to *in vivo* studies in order to assess the feasibility of further, possibly expensive, experimentation.

## TRANSECT Web-application
A web-accessible version is freely available for use [here](http://203.101.229.190/analysis/home)

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

NOTE: TRANSECT requires and depends on numerous packages and applications. For this reason it is recommended for most cases, to install TRANSECT within a Conda environment. This is not only easier and cleaner but also saves a significant amount of time (*from approx. 60 minutes to less than 10*). For those not wanting to use Conda, please see the instructions for an alternate method in the TRANSECT manual.

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

There are many wikis on how to install Conda on Ubuntu, [here](https://docs.anaconda.com/miniconda/) is just one.

*(approx. 1min)*

```sh
### follow the instructions outlined in the link above which should look something like this
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh

### don't forget to initialize the bash shell
~/miniconda3/bin/conda init bash
```

Create the TRANSECT Conda environment. In this step, we will create an environment with all the tools and dependencies required to run TRANSECT. 

*(approx. 5-10mins)*

```sh
### change into the top directory of the downloaded folder (TRANSECT) and navigate to the INSTALL folder
cd <path to>/TRANSECT/INSTALL

### First, run the conda install script 
./TRANSECT_conda_install.sh

### Next, upon succesful completion of the previous step, activate the newly created environment
conda activate TRANSECT

### Finally, once in the TRANSECT environment, run the post installation script to complete the setup
./TRANSECT_post_conda_install.sh

# Note: you should reactivate the TRANSECT environment at this point 
conda deactivate
conda activate TRANSECT
```

And that's it! You should now have all the necessary applications and dependencies to run TRANSECT. Please note, just like any virtual environment you are required to activate the TRANSECT environment in order to use the application. You can deactivate at will when not in use. 

## Usage

<span id="#usage"></span> 

TRANSECT has two main operations; **Prepare** and **Analyse**. 

**Prepare** is a process that retrieves the raw data from online repositories and prepares it (if required) for analysis. TRANSECT comes bundled with three different prepare scripts, one each for RECOUNT3, GDC-TCGA and GTEx data. 

Example prepare commands;

```sh
### NOTE: use the -h parameter with any of the following scripts to see the full help menu. ie.
# bin/GDC_TCGA_prepare_directories.sh -h

### ALSO NOTE: you don't necessarily need all of these, all at once. If your just tying this application start with RECOUNT3 (recommended) and skip the others for now.

### RECOUNT3 (approx. 5mins)
# change into the top directory of TRANSECT
cd <path to>/TRANSECT
# run the RECOUNT3 prepare script for TCGA-BRCA
bin/R3_prepare_directories.sh -p BRCA

### GDC-TCGA
# change into the top directory of TRANSECT
cd <path to>/TRANSECT
# run the GDC-TCGA prepare script for TCGA-BRCA requiring here, a complete download (mRNA and miRNA)
bin/GDC_TCGA_prepare_directories.sh -p TCGA-BRCA -a

### GTEx
# change into the top directory of TRANSECT
cd <path to>/TRANSECT
# run the GTEx prepare script (NOTE: GTEx data is retreived not by tissue as above, but as a single package for all tissues)
bin/GTEx_prepare_directories.sh -a
```

Be aware that some of these collections are large and require substantial disk space. They will take a considerable amount of time to download and process too. For example, downloading and processing TCGA-BRCA takes just over 30 minutes (using a high speed network connection and an up to date workstation) and requires more than 14GB of disk space (most of which can and by default is, deleted afterwards). In comparison, TCGA-LAML takes less than 5 minutes to retrieve and less than 2GB of disc space. Preparation for both GDC-TCGA and the RECOUNT3 data is done individually by tissue type but can also be done in batch mode (see the relevant script help menu for instructions). Preparing GTEx data on the other hand, retrieves in bulk all tissue types in a single table before separating them into individual files based on tissue type (again, see the help menu for more details). All downloaded data is stored in "TRANSECT/data/<GDC|GTEx|RECOUNT3>" in appropriate folders.

**Analyse** is a process that uses the prepared public data from above conducts the stratified differential expression and produces all the outputs. Like with the prepare operation, TRANSECT comes bundled with three analyse scripts, one each for RECOUNT3, GDC-TCGA and GTEx data.

Example analyse commands;

```sh
### NOTE: use the -h parameter with any of the following scripts to see the full help menu. ie.
# bin/GDC_TCGA_analyse_GOI.sh -h

### ALSO NOTE: You need to have the appropriate DB installed (previous step) for the following commands to work. If you only installed the RECOUNT3 BRCA DB, only run the RECOUNT3 BRCA test.

### RECOUNT3 (approx. 5-10mins)
# change into the top directory of TRANSECT
# TRANSECT saves output in the run folder so best to create one specifically for each run 
cd <path to>/TRANSECT/output/RECOUNT3
mkdir -p ESR1_BRCA_test
cd ESR1_BRCA_test
# Now, run the RECOUNT3 prepare script on the BRCA data, investigating the gene ESR1, with all outputs.
<path to>/TRANSECT/bin/R3_analyse_GOI.sh -p BRCA -g ESR1 -s mRNA -t 3 -a

### GDC-TCGA
# TRANSECT saves output in the run folder so best to create one specifically for each run 
cd <path to>/TRANSECT/output/GDC
mkdir -p ESR1_BRCA_test
cd ESR1_BRCA_test
# Now, run the GDC-TCGA prepare script on the BRCA data, investigating the gene ESR1, with all outputs
<path to>/TRANSECT/bin/GDC_TCGA_analyse_GOI.sh -p TCGA-BRCA -g ESR1 -s mRNA -t 3 -a

### GTEx
# TRANSECT saves output in the run folder so best to create one specifically for each run 
cd <path to>/TRANSECT/output/GTEx
mkdir -p ESR1_Breast_test
cd ESR1_Breast_test
# Now, run the GTEx prepare script on the Breast data, investigating the gene ESR1, with all outputs
<path to>/TRANSECT/bin/GTEx_analyse_GOI.sh -p BREAST -g ESR1 -s mRNA -t 5 -a
```

## Output

<span id="#output"></span>

...

## Manual

<span id="#manual"></span>

...

## Contributors

...

## Licence

TRANSECT is made available under the terms of the MIT license, a copy of which is included in the distribution in the LICENSE file.

GSEA is made available under the terms of a BSD-style license, a copy of which is included in the GSEA folder in the LICENSE.txt file. See that file for exact terms and conditions.

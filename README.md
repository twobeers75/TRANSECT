<img title="TRANSECT" alt="TRANSECT" src="https://github.com/twobeers75/TRANSECT/blob/main/TRANSECT_t.png">

#
| [Installation](#installation) | [Usage](#usage) | [Output](#output) | [Manual](#manual) |


## TRANSECT Introduction

TRANSECT is an application designed solely to run in a Linux Terminal. In brief, TRANSECT stratifies large cohort data into user defined subsets which are subsequently compared one to the other. This version of TRANSECT is tuned specifically for the interrogation of gene/s from publicly available tissue-specific gene expression data.

Currently, there are many global efforts to collect and collate tissue-specific gene expression (RNA-sequencing) data from large cohorts of diseased and non-diseased subjects for public use (GTEx, TCGA, …, to name only a few). These data sets are used to bolster and aid in the investigation of gene specific expression and regulation as well as several other biological queries. One such application that is currently underutilised and (to our knowledge) not currently offered as an accessible public resource, is the stratification and subsequent differential expression of single or composite genes in order to assess "transcriptome state". This type of analysis has only in recent times become feasible thanks to the large number of participants in cohort studies and, will only intensify in power and usefulness as these collections continue to grow in number and diversity. The results from this type of analysis can be used to rapidly investigate transcriptome state within natural physiological expression levels without economic burden. This can be applied before mass-parallel laboratory investigations and/or between the transition from *in vitro* to *in vivo* studies in order to assess the feasibility of further, possibly expensive, experimentation.

## TRANSECT Web-application
A web-accessible version is freely available for use [here](https://transect.au)

## Version

24.03

## Authors

....

## Hardware requirements


64 bit Linux (developed solely on Ubuntu 22.04)

10GB RAM (12+ recommended)

## Software requirements

- Python 3
  - see pip_requirements.txt for a complete list of required python packages
- R
  - see r_requirements.txt for a complete list of required R packages
- Java
  - preinstalled with Ubuntu 22.04 (usually, we'll make sure during this installation)


## Installation
<span id="#installation"></span>

NOTE: TRANSECT requires and depends on numerous packages and applications. These take some time to install if not already present. A fresh install on a vanilla Ubuntu 22.04 can take 30-45mins depending on the PC and network speeds. 

To start, clone the repo

```sh
git clone https://github.com/twobeers75/TRANSECT.git
```

Or, download the package and extract it

Browse to https://github.com/twobeers75/TRANSECT, click on the green "Code" button followed by "Download ZIP" (note the download location). 
Find the downloaded ZIP file and move it if required before extracting the contents and renaming the folder

```sh
unzip TRANSECT-main.zip
mv TRANSECT-main TRANSECT
```

Install python3 pip, java if required, and other TRANSECT dependencies (NOTE: Python3 comes preinstalled with Ubuntu)

*(approx. 1-2min)*
```sh
### change into the top directory of the downloaded folder (TRANSECT)
cd <path to>/TRANSECT

### install pip and other deb requirements
sudo apt install python3-pip default-jre libfontconfig1-dev libcurl4-openssl-dev libssl-dev libxml2-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev pandoc

### install python modules
python3 -m pip install -r pip_requirements.txt
```

Install R, the "pacman" package and Bioconductor specific packages. You can skip this step if you already have R.

There are many wikis on how to install R on Ubuntu, [here](https://www.digitalocean.com/community/tutorials/how-to-install-r-on-ubuntu-22-04) is just one (specifically for Ubuntu 22.04)

*(approx. 1min)*
```sh
# follow the instructions outlined in the link above which should look something like this
wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo gpg --dearmor -o /usr/share/keyrings/r-project.gpg
echo "deb [signed-by=/usr/share/keyrings/r-project.gpg] https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" | sudo tee -a /etc/apt/sources.list.d/r-project.list
sudo apt update
sudo apt install r-base
```

Start R from the terminal and install pacman and devtools. Follow the prompts and choose (if asked) to install these packages into a personal library.

Once you enter the R shell you should see printed out in the terminal a number of lines about the R version and licences followed by a ">" symbol. I have used this symbol below to indicate that you need to be in the R shell to run these commands but, you can't copy the ">" symbol too. It won't work. 

*(approx. 25mins)*
```
### start R
R
> install.packages(c("pacman","devtools"))
# Note: maybe wise here to go get a coffee as the previous command takes quite some time to finish! (approx. 15mins) 

### whilst still in the R environment, load devtools and install rlogging
> library("devtools")
> install_github("https://github.com/mjkallen/rlogging.git")

### also install required Bioconductor packages
> if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
> BiocManager::install(version = "3.19")
> BiocManager::install(c("edgeR","Glimma","DEFormats"))
# Note: probably time for another coffee. Sorry! (approx. 10mins)

### once successfully completed you can quit R, no need to save the workspace.
# No more coffee for you today ;-)
>  q()
```



NOTE: TRANSECT requires many additional R packages however these are all installed on demand the first time (and only the first time) you run each one of the different TRANSECT commands. Please keep this in mind on your first run as it will take substantially longer compared to all subsequent runs.

## Usage

<span id="#usage"></span> 

TRANSECT has two main operations; **Prepare** and **Analyse**. 

**Prepare** is a process that retrieves the raw data from online repositories and prepares it (if required) for analysis. TRANSECT comes bundled with three different prepare scripts, one each for RECOUNT3, GDC-TCGA and GTEx data. 

Example prepare commands;

```sh
### NOTE: use the -h parameter with any of the following scripts to see the full help menu. ie.
# bin/GDC_TCGA_prepare_directories.sh -h

### ALSO NOTE: you don't necessarily need all of these, all at once. If your just tying this application start with RECOUNT3 (recommended) and skip the others for now.

### RECOUNT3 (approx. 10mins on first run)
# change into the top directory of TRANSECT
cd <path to>/TRANSECT
# run the RECOUNT3 prepare script for BRCA requiring a complete download. If this is a first time download, the script will autmatically install all requirements so it will take a while but, next time will be quick! 
bin/R3_prepare_directories.sh -p BRCA

### GDC-TCGA
# change into the top directory of TRANSECT
cd <path to>/TRANSECT
# run the GDC-TCGA prepare script for BRCA requiring a complete download
bin/GDC_TCGA_prepare_directories.sh -p TCGA-BRCA -a

### GTEx
# change into the top directory of TRANSECT
cd <path to>/TRANSECT
# run the GTEx prepare script requiring a complete download
bin/GTEx_prepare_directories.sh -a
```

Be aware that some of these collections are large and require substantial disk space. They will take a considerable amount of time to download and process too. For example, downloading and processing TCGA-BRCA takes just over 30 minutes (using a high speed network connection and an up to date workstation) and requires more than 14GB of disk space (most of which can be deleted afterwards). In comparison, TCGA-LAML takes less than 5 minutes to retrieve and less than 2GB of disc space. Preparation for both GDC-TCGA and the RECOUNT3 data is done individually by tissue type but can also be done in batch mode (see the relevant script help menu for instructions). Preparing GTEx data on the other hand, retrieves in bulk all tissue types in a single table before separating them into individual files based on tissue type (again, see the help menu for more details). All downloaded data is stored in "TRANSECT/data/<GDC|GTEx|RECOUNT3>" in appropriate folders.

**Analyse** is a process that uses the prepared public data from above conducts the stratified differential expression and produces all the outputs. Like with the prepare operation, TRANSECT comes bundled with three analyse scripts, one each for RECOUNT3, GDC-TCGA and GTEx data.

Example analyse commands;

```sh
### NOTE: use the -h parameter with any of the following scripts to see the full help menu. ie.
# bin/GDC_TCGA_analyse_GOI.sh -h

### ALSO NOTE: You need to have the appropriate DB installed (previous step) for the following commands to work. If you only installed the RECOUNT3 BRCA DB, only run the RECOUNT3 BRCA test.

### RECOUNT3 (approx. 10mins on first run)
# change into the top directory of TRANSECT
# TRANSECT saves output in the run folder so best to create one specifically for each run 
cd <path to>/TRANSECT/output/RECOUNT3
mkdir -p ESR1_BRCA_test
cd ESR1_BRCA_test
# Now, run the RECOUNT3 prepare script on the BRCA data, investigating the gene ESR1, with all outputs. Again, if this is a first time analysis, the script will autmatically install all requirements so it will take a while but, next time will be quick!
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
<path to>/TRANSECT/bin/GTEx_analyse_GOI.sh -p Breast -g ESR1 -s mRNA -t 5 -a
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

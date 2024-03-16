# SCA - Stratified Cohort Analysis
SCA is an application designed solely to run in a Linux Terminal. In brief, SCA stratifies large cohort data into user defined subsets which are subsequently compared one to the other. This version of SCA is tuned specifically for the interrogation of gene/s from publicly available tissue-specific gene expression data.

Currently, there are many global efforts to collect and collate tissue-specific gene expression (RNA-sequencing) data from large cohorts of diseased and non-diseased subjects for public use (GTEx, TCGA, …, to name only a few). These data sets are used to bolster and aid in the investigation of gene specific expression and regulation as well as several other biological queries. One such application that is currently underutilised and (to our knowledge) not currently offered as an accessible public resource, is the stratification and subsequent differential expression of single or composite genes in order to assess "transcriptome fate". This type of analysis has only in recent times become feasible thanks to the large number of participants in cohort studies and, will only intensify in power and usefulness as these collections continue to grow in number and diversity. The results from this type of analysis can be used to rapidly investigate transcriptome state within natural physiological expression levels without economical burden. This can be applied before mass-parallel laboratory investigations and/or between the transition from *in vitro* to *in vivo* studies in order to assess the feasibility of further, possibly expensive, experimentation.

## Version

24.03

## Authors

....

## Hardware requirements

64 bit Linux (tested solely on Ubuntu 22.04)

>10GB RAM (12+ recommended)

## Software dependencies

- Python 3
  - see pip_requirements.txt for a complete list of required python packages
- R
  - see r_requirements.txt for a complete list of required R packages

## Installation
NOTE: SCA requires and depends on numerous packages and applications. These take some time to install if not already present. A fresh install on a vanilla Ubuntu 22.04 can take 30-45mins depending on the PC and network speeds. 

Clone the repo

```sh
git clone https://github.com/twobeers75/SCA.git
```

Or, download the package and extract it
Browse to https://github.com/twobeers75/SCA, click on the green "Code" button followed by "Download ZIP" (note the download location). 
Find the downloaded ZIP file and move it if required before extracting the contents and renaming the folder
```sh
unzip SCA-main.zip
mv SCA-main SCA
```

Install python3 pip and all python dependencies (NOTE: Python 3 comes preinstalled with Ubuntu)

```sh
### change into the top directory of the downloaded folder (SCA)
cd <path to>/SCA

### install pip and other requirements
sudo apt install python3-pip libfontconfig1-dev libcurl4-openssl-dev libssl-dev libxml2-dev

### install python modules
python3 -m pip install -r pip_requirements.txt
```

Install R, the "pacman" package and Bioconductor specific packages.

```sh
### There are many wikis on how to install R on Ubuntu, here is just one
### https://www.digitalocean.com/community/tutorials/how-to-install-r-on-ubuntu-22-04
# follow the instructions outlined in the link above which should look something like this
wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo gpg --dearmor -o /usr/share/keyrings/r-project.gpg
echo "deb [signed-by=/usr/share/keyrings/r-project.gpg] https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" | sudo tee -a /etc/apt/sources.list.d/r-project.list
sudo apt update
sudo apt install r-base

### Start R from the terminal and install pacman and devtools (follow the prompts and choose (if asked) to install into a personal library)
R
> install.packages(c("pacman","devtools"))

### Load devtools and install rlogging
> library("devtools")
> install_github("https://github.com/mjkallen/rlogging.git")

### Also install required Bioconductor packages (these take quite some time to install!)
> if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
> BiocManager::install(version = "3.18")
> BiocManager::install(c("edgeR","Glimma","DEFormats"))

# Once successfully completed you can quit R
>  q()
```



NOTE: SCA requires many additional R packages however these are all installed on demand the first time (and only the first time) you run each one of the different SCA commands. Please keep this in mind on your first run as it will take substantially longer compared to all subsequent runs.

## Usage

 SCA has two main operations; **Prepare** and **Analyse**. 

**Prepare** is a process that retrieves the raw data from online repositories and prepares it (if required) for analysis. SCA comes bundled with three different prepare scripts, one each for GDC-TCGA, GTEx and RECOUNT3 data. 

Example prepare commands;

```sh
### NOTE: use the -h parameter with any of the following scripts to see the full help menu. ie.
# bin/GDC_TCGA_prepare_directories.sh -h

### RECOUNT3
# change into the top directory of SCA
cd <path to>/SCA
# run the RECOUNT3 prepare script for BRCA requiring a complete download
bin/R3_prepare_directories.sh -p BRCA

### GDC-TCGA
# change into the top directory of SCA
cd <path to>/SCA
# run the GDC-TCGA prepare script for BRCA requiring a complete download
bin/GDC_TCGA_prepare_directories.sh -p TCGA-BRCA -a

### GTEx
# change into the top directory of SCA
cd <path to>/SCA
# run the GTEx prepare script requiring a complete download
bin/GTEx_prepare_directories.sh -a
```

Be aware that some of these collections are large and require substantial disk space. They will take a considerable amount of time to download and process too. For example, downloading and processing TCGA-BRCA takes just over 30 minutes (using a high speed network connection and an up to date workstation) and requires more than 14GB of disk space (most of which can be deleted afterwards). In comparison, TCGA-LAML takes less than 5 minutes to retrieve and less than 2GB of disc space. Preparation for both GDC-TCGA and the RECOUNT3 data is done individually by tissue type but can also be done in batch mode (see the relevant script help menu for instructions). Preparing GTEx data on the other hand, retrieves in bulk all tissue types in a single table before separating them into individual files based on tissue type (again, see the help menu for more details). All downloaded data is stored in "SCA/data/<GDC|GTEx|RECOUNT3>" in appropriate folders.

**Analyse** is a process that uses the prepared public data from above and produces stratified differential expression and other associated outputs. Like with the prepare operation, SCA comes bundled with three analyse scripts, one each for GDC-TCGA, GTEx and RECOUNT3 data.

Example analyse commands;

```sh
### NOTE: use the -h parameter with any of the following scripts to see the full help menu. ie.
# bin/GDC_TCGA_analyse_GOI.sh -h

### RECOUNT3
# change into the top directory of SCA
# SCA saves output in the run folder so best to create one specifically for each run 
cd <path to>/SCA/output/RECOUNT3
mkdir -p ESR1_BRCA_test
cd ESR1_BRCA_test
# Now, run the RECOUNT3 prepare script on the BRCA data, investigating the gene ESR1, with all outputs
<path to>/SCA/bin/R3_analyse_GOI.sh -p BRCA -g ESR1 -s mRNA -t 3 -a

### GDC-TCGA
# SCA saves output in the run folder so best to create one specifically for each run 
cd <path to>/SCA/output/GDC
mkdir -p ESR1_BRCA_test
cd ESR1_BRCA_test
# Now, run the GDC-TCGA prepare script on the BRCA data, investigating the gene ESR1, with all outputs
<path to>/SCA/bin/GDC_TCGA_analyse_GOI.sh -p TCGA-BRCA -g ESR1 -s mRNA -t 3 -a

### GTEx
# SCA saves output in the run folder so best to create one specifically for each run 
cd <path to>/SCA/output/GTEx
mkdir -p ESR1_Breast_test
cd ESR1_Breast_test
# Now, run the GTEx prepare script on the Breast data, investigating the gene ESR1, with all outputs
<path to>/SCA/bin/GTEx_analyse_GOI.sh -p Breast -g ESR1 -s mRNA -t 5 -a
```

## Output

...

## Manual

...

## Contributors

...

## Licence

...


#!/usr/bin/env bash

##########################################################################################################
#### INSTALL TRANSECT ENVIRONMENT CONDA
#### JToubia - January 2024
#########

### Check conda is installed on system
if [[ $(command -v conda) ]]; then 
	echo "Conda installed"
	conda_year=`conda --version | cut -f 2 -d " " | cut -f 1 -d "."`
	if [ "$conda_year" -lt 24 ]; then
		echo "Conda version is too old, requires conda > 24.1.0"
		echo "Please use the following command to update to latest release"
		echo "conda update conda"
		exit 1
	else
		echo "Conda version is compatible"
	fi
else
	echo "It appears Conda is not installed but is required"
	echo "Before proceeding you will first need to install conda"
	echo "See the link below for instructions on how to install Miniconda for Linux" 
	echo "https://docs.anaconda.com/miniconda/"
	echo ""
	echo "Afterwards, please return to this step and try again"
	echo "NOTE: Don't forget to initialize your newly-installed Miniconda in bash"
	echo "(~/miniconda3/bin/conda init bash)" 
	echo ""
	echo "and restart your terminal too!"
	echo ""
	exit 1
fi

### Create TRANSECT environment
echo "Creating TRANSECT conda environment"
echo ""
conda create --yes -n TRANSECT -c bioconda -c conda-forge -c anaconda \
		fontconfig \
		libcurl \
		openssl \
		libxml2 \
		harfbuzz \
		fribidi \
		freetype \
		libpng \
		libtiff \
		libjpeg-turbo \
		pandoc \
		git \
		zip \
	python \
		pip \
		glob2 \
		numpy \
		pandas \
		pyarrow \
		requests \
	r-base \
		r-pacman \
		r-devtools \
		r-biocmanager \
		r-data.table \
		r-dplyr \
		r-tidyr \
		r-tibble \
		r-ggplot2 \
		r-ggforce \
		r-gplots \
		r-calibrate \
		r-rcolorbrewer \
		r-plotly \
		r-htmlwidgets \
		r-webgestaltr \
		r-zip \
		bioconductor-edger \
		bioconductor-glimma \
		bioconductor-deformats \
		bioconductor-recount3=1.12.0 \
		bioconductor-deformats \
	anaconda::openjdk



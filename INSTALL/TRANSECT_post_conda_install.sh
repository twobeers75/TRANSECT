#!/usr/bin/env bash

##########################################################################################################
#### INSTALL TRANSECT REQUIRED rlogging and set R_ZIPCMD env variable
#### JToubia - January 2024
#########
### Check to see we are in the TRANSECT environment
current_conda_env="${CONDA_DEFAULT_ENV}"
if [[ ${current_conda_env} != "TRANSECT" ]]; then
	echo "Appears as if you are not in the TRANSECT conda environment"
	echo "Run the following command before trying this script again"
	echo ""
	echo "conda activate TRANSECT"
	echo ""
	exit 1
fi

### Check we are in the INSTALL folder and record bin dir location
workdir=$(pwd)
bindir=$(dirname "$workdir")/bin
if [ ! -d "${bindir}" ]; then
	echo "Appears that we can't find the TRANSECT bin directory where it should be"
	echo "If the TRANSECT INSTALL folder has been moved try moving it back to its original location"
	echo "before rerunning this script"
	exit 1
fi

### Install rlogging
echo "Installing required R package - rlogging"
## get rlogging
git clone https://github.com/mjkallen/rlogging.git

## install in R
R CMD INSTALL rlogging

## remove rlogging
rm -rf rlogging
echo ""

### Issues with R accessing zip in Conda, need to check if environment variable for zip is set (TODO: should actually do this all in a .yml setup file)
echo "Checking that the R_ZIPCMD environment variable is set in this Conda environment"
zip_exe=$(which zip)
r_zip_exe=$(Rscript -e 'Sys.getenv("R_ZIPCMD","zip")')
if [[ ${r_zip_exe} != *\/* ]]; then
	echo "Does not appear to be set, setting it now"
	conda env config vars set R_ZIPCMD="${zip_exe}"
fi

### Finally, add TRANSECT bin directory to PATH in conda environment
## need to save these scripts in the right place using $CONDA_PREFIX
## add activate script
echo '#!/usr/bin/env bash' > ${CONDA_PREFIX}/etc/conda/activate.d/transect_activate.sh
echo "bindir=${bindir}" >> ${CONDA_PREFIX}/etc/conda/activate.d/transect_activate.sh
echo 'export PATH="$bindir:$PATH"' >> ${CONDA_PREFIX}/etc/conda/activate.d/transect_activate.sh
## add deactivate script
echo '#!/usr/bin/env bash' > ${CONDA_PREFIX}/etc/conda/deactivate.d/transect_deactivate.sh
echo "bindir=${bindir}" >> ${CONDA_PREFIX}/etc/conda/deactivate.d/transect_deactivate.sh
echo 'export PATH=$(echo $PATH | tr ":" "\n" | grep -v $bindir | xargs | tr " " ":")' >> ${CONDA_PREFIX}/etc/conda/deactivate.d/transect_deactivate.sh

### Done sign off with message
echo ""
echo "Environment setup complete!"
echo "You are now required to reactivate this environment before use;"
echo ""
echo "conda deactivate"
echo "conda activate TRANSECT"
echo ""


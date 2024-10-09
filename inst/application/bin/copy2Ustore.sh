#!/bin/bash

# Set default values for filename
filename="summarizedExperiment.rds"

# Check if arguments were passed, if so, override defaults
project=${1}   # First argument for project, default to 33485
workunit=${2}  # Second argument for workunit, must be the folder containing the rds file.
if [ ! -z "$3" ]; then
    filename=$3  # Optional third argument for filename
fi

year=$(date +"%Y")
yearmonth=$(date +"%Y-%m")
yearmonthday=$(date +"%Y-%m-%d")

# Create paths similar to standard and stage the data on your scratch on the 35
url="uStore/p${project}/bfabric/Proteomics/SummarizedExperiment/${year}/${yearmonth}/${yearmonthday}/workunit_${workunit}/$filename"
echo "$url"

tmppath="wolski@fgcz-r-035:/srv/www/htdocs/public/$url"

# Remove hostname from path and get directory name
path_without_hostname=$(echo "$tmppath" | sed 's/^[^:]*://')
path_without_filename=$(dirname "$path_without_hostname")
echo "$path_without_filename"

# Find the first .rds file in the workunit directory
RDSFILE=$(find "$workunit" -type f -name "*.rds" | head -n 1)
echo "$RDSFILE"
ls -l "$RDSFILE"

# Create directory on remote server
ssh wolski@fgcz-r-035 "mkdir -p $path_without_filename"

# Copy the .rds file to the remote path
scp -r "$RDSFILE" "$tmppath"

# Generate URL and append it to upload.txt
URL0="https://fgcz-shiny.uzh.ch/exploreDE/?data=../uStore/"
OUTPUTURL=$url
URL1=$(echo "$OUTPUTURL" | sed 's/^.*\/\(p[0-9][0-9]*\/.*rds\)$/\1/') || { echo "could not derive valid url link argument"; exit 1; }


echo "${URL0}${URL1}"
echo "${URL0}${URL1}" >> upload.txt

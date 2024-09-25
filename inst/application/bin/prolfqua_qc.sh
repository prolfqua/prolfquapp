#!/bin/bash


# Function to extract the directory of the dataset argument
get_dataset_dir() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            -d|--dataset)
                dataset_path="$2"
                dataset_dir=$(dirname "$dataset_path")
                return 0
                ;;
        esac
        shift
    done
    return 1  # Return error if no dataset argument is found
}


# Get the path to the installed R package
PACKAGE_PATH=$(Rscript --vanilla -e "cat(system.file(package = 'prolfquapp'))")

# Build the path to the R script
R_SCRIPT_PATH="${PACKAGE_PATH}/application/CMD_QUANT_QC.R"


# Default log file path
log_file="prolfqua_logMemUsage_qc.log"

# Extract dataset directory if -d or --dataset argument is provided
if get_dataset_dir "$@"; then
    log_file="${dataset_dir}/prolfqua_logMemUsage_dataset.log"
    echo $log_file
    # Create the directory if it does not exist
    mkdir -p "$dataset_dir"
fi

# Check if the R script exists
if [[ -f "$R_SCRIPT_PATH" ]]; then
    echo "Rscript --vanilla \"$R_SCRIPT_PATH\" \"$@\""
    Rscript --vanilla "$R_SCRIPT_PATH" "$@" &
    PID=$!
    while ps -p $PID > /dev/null
    do
      ps -p $PID -o %cpu,%mem,etime,rss,utime,stime,cputime >> "$log_file"
      sleep 1
    done
else
    echo "Error: R script not found at '$R_SCRIPT_PATH'"
    exit 1
fi

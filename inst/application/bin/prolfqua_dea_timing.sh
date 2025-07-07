#!/bin/bash


# Function to extract the directory of the dataset argument
get_indir() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            -i|--indir)
                indir_dir="$2"
                return 0
                ;;
        esac
        shift
    done
    return 1  # Return error if no dataset argument is found
}


# Function to rename the log file if it exists by adding a number (1, 2, 3, ...)
rename_log_file_if_exists() {
    local file="$1"
    local counter=1
    local new_file="${file%.log}_$counter.log"

    while [[ -f "$new_file" ]]; do
        counter=$((counter + 1))
        new_file="${file%.log}_$counter.log"
    done

    if [[ -f "$file" ]]; then
        mv "$file" "$new_file"
        echo "Renamed existing log file to $new_file"
    fi
}




# Default log file path
log_file="prolfqua_logMemUsage_dea.log"

# Extract dataset directory if -d or --dataset argument is provided
if get_indir "$@"; then
    log_file="${indir_dir}/prolfqua_logMemUsage_dea.log"
    echo $log_file
    # Create the directory if it does not exist
    mkdir -p "$indir_dir"
fi

rename_log_file_if_exists $log_file


# Get the path to the installed R package
PACKAGE_PATH=$(Rscript --vanilla -e "cat(system.file(package = 'prolfquapp'))")
# Build the path to the R script
R_SCRIPT_PATH="${PACKAGE_PATH}/application/CMD_DEA.R"


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

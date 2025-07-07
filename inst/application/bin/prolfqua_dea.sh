#!/bin/bash

set -euxo pipefail

# Get the path to the installed R package
PACKAGE_PATH=$(Rscript --vanilla -e "cat(system.file(package = 'prolfquapp'))")
# Build the path to the R script
R_SCRIPT_PATH="${PACKAGE_PATH}/application/CMD_DEA.R"


if [[ -f "$R_SCRIPT_PATH" ]]; then
    echo "Rscript --vanilla \"$R_SCRIPT_PATH\" \"$@\""
    Rscript --vanilla "$R_SCRIPT_PATH" "$@"
else
    echo "Error: R script not found at '$R_SCRIPT_PATH'"
    exit 1
fi

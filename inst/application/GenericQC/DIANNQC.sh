#!/bin/bash

# Function to display help
function show_help {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -h, --help                 Show this help message and exit"
    echo "  -i STRING, --indir=STRING         folder containing required files"
    echo "  -p STRING, --projectId=STRING     your project identifier"
    echo "  -w STRING, --workunitId=STRING    workunit identifier"
    echo "  -d STRING, --dataset=STRING       name of annotation"
    echo "  -o STRING, --output=STRING        folder to write the results to"
    echo "  --libPath=STRING          (optional) R library path"
}

# Default parameters
indir=""
projectId=""
workunitId=""
dataset=""
output=""
libPath=""

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -h|--help) show_help; exit 0 ;;
        -i|--indir) indir="$2"; shift ;;
        -p|--projectId) projectId="$2"; shift ;;
        -w|--workunitId) workunitId="$2"; shift ;;
        -d|--dataset) dataset="$2"; shift ;;
        -o|--output) output="$2"; shift ;;
        --libPath) libPath="$2"; shift ;;
        --) shift; break ;;
        -*) echo "Error: Unsupported flag $1" >&2; show_help; exit 1 ;;
        *) break ;;
    esac
    shift
done

# Run the R script with the provided arguments
Rscript --vanilla DIANN_QC_V2.R \
    --indir="$indir" \
    --projectId="$projectId" \
    --workunitId="$workunitId" \
    --dataset="$dataset" \
    --output="$output" \
    $( [[ -n $libPath ]] && echo "--libPath=$libPath" )

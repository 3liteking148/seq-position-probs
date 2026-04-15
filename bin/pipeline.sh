#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <hmm_file> <fa_file>"
    exit 1
fi

HMM_FILE="$1"
FA_FILE="$2"

# Create a secure temporary file for the FASTA output
TEMP_FA=$(mktemp /tmp/temp_stderr_XXXXXX.fa)

# Set a trap to ensure the temporary file is deleted when the script exits
# (whether it finishes successfully or is interrupted)
trap 'rm -f "$TEMP_FA"' EXIT

# Run bathsearch and redirect only stderr (2>) to the temporary file.
# Standard output (stdout) will still print to the terminal unless redirected elsewhere.
$SCRIPT_DIR/bathsearch --cpu 1 "${HMM_FILE}.filter.bhmm" "$FA_FILE" 2> "$TEMP_FA"

# Check if the filter command succeeded (optional, but good practice)
if [ $? -ne 0 ]; then
    echo "Error: bathsearch encountered an issue."
    exit 1
fi

# Run dummer using the original HMM file and the captured temporary FASTA file
$SCRIPT_DIR/dummer -v "$HMM_FILE" "$TEMP_FA"
cp $TEMP_FA /mnt/tmp/test.txt
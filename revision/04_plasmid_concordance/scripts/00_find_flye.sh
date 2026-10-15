#!/bin/bash
# Helper: list candidate Flye assembly FASTAs for one isolate, to set FLYE_GLOB in config.sh.
# Usage: bash scripts/00_find_flye.sh 6817
source "$(dirname "$0")/../config.sh"
ID=${1:?give an isolate number}
find "$NEWCRE" -path "*${ID}*" \( -name "*.fasta" -o -name "*.fa" \) 2>/dev/null | grep -i -E "flye|medaka|consensus|assembly" | sort

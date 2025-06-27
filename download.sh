#!/bin/bash

set -euo pipefail

MANIFEST_FILE="neale_lab_manifest_20180731.tsv"

show_help() {
    cat << EOF
Usage: $0 [--phenotype <description> | --phenotype_code <code>] --output <output_directory>

Download GWAS data for a specific phenotype from the Neale lab UK Biobank dataset.
Only downloads both_sexes and irnt normalized files.

Options:
  --phenotype       Phenotype description (e.g., "Calcium", "Food weight")
  --phenotype_code  Phenotype code (e.g., 100024_irnt, 100001_irnt)
  --output          Output directory name (will be created if it doesn't exist)
  --help            Show this help message

Examples:
  $0 --phenotype "Calcium" --output calcium
  $0 --phenotype_code 100024_irnt --output calcium
  $0 --phenotype "Food weight" --output food_weight
EOF
}

# Parse command line arguments
PHENOTYPE=""
PHENOTYPE_CODE=""
OUTPUT_DIR=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --phenotype)
            PHENOTYPE="$2"
            shift 2
            ;;
        --phenotype_code)
            PHENOTYPE_CODE="$2"
            shift 2
            ;;
        --output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --help)
            show_help
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            show_help
            exit 1
            ;;
    esac
done

# Validate required arguments
if [[ -z "$PHENOTYPE" && -z "$PHENOTYPE_CODE" ]]; then
    echo "Error: Either --phenotype or --phenotype_code is required" >&2
    show_help
    exit 1
fi

if [[ -n "$PHENOTYPE" && -n "$PHENOTYPE_CODE" ]]; then
    echo "Error: Cannot specify both --phenotype and --phenotype_code" >&2
    show_help
    exit 1
fi

if [[ -z "$OUTPUT_DIR" ]]; then
    echo "Error: --output is required" >&2
    show_help
    exit 1
fi

# Check if manifest file exists
if [[ ! -f "$MANIFEST_FILE" ]]; then
    echo "Error: Manifest file $MANIFEST_FILE not found" >&2
    exit 1
fi

# Search for the phenotype in the manifest (both_sexes and _irnt only)
if [[ -n "$PHENOTYPE_CODE" ]]; then
    # Search by phenotype code (column 1)
    WGET_CMD=$(awk -F'\t' -v phenotype_code="$PHENOTYPE_CODE" '
        NR > 1 && $1 == phenotype_code && $4 == "both_sexes" && $1 ~ /_irnt$/ {
            print $6
            exit
        }
    ' "$MANIFEST_FILE")
    SEARCH_TERM="$PHENOTYPE_CODE"
else
    # Search by phenotype description (column 2)
    WGET_CMD=$(awk -F'\t' -v phenotype_desc="$PHENOTYPE" '
        NR > 1 && $2 == phenotype_desc && $4 == "both_sexes" && $1 ~ /_irnt$/ {
            print $6
            exit
        }
    ' "$MANIFEST_FILE")
    SEARCH_TERM="$PHENOTYPE"
fi

if [[ -z "$WGET_CMD" ]]; then
    echo "Error: Phenotype '$SEARCH_TERM' not found in manifest or does not match criteria (both_sexes + _irnt)" >&2
    echo "Available _irnt phenotypes for both_sexes (Random selection of 10 rows):"
    awk -F'\t' 'NR > 1 && $4 == "both_sexes" && $1 ~ /_irnt$/ { print "  " $1 " - " $2 }' "$MANIFEST_FILE" | shuf | head -10
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

echo "Downloading $PHENOTYPE to directory: $OUTPUT_DIR"
echo "Command: $WGET_CMD"

# Execute the wget command
eval "$WGET_CMD"

echo "Download completed successfully!"
echo "Files in $OUTPUT_DIR:"
ls -la

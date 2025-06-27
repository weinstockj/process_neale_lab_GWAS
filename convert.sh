#!/bin/bash

input="$1"
output="$2"

if [ -z "$input" ]; then
    echo "Usage: $0 <input_file> <output_file>" >&2
    exit 1
fi

# dbSNP directory path from project documentation
dbsnp_dir="/home/jwein22/weinstocklab/projects/resources/dbSNP_to_parquet"
chain="hg19ToHg38.over.chain"
VARIANT_FILE="/beegfs/labs/weinstocklab/projects/PRS/LD_REFERENCE_PANEL/output/tsv/variant_list_EUR.tsv"

echo "Converting GWAS data using convert.py..."

pixi run python convert.py -s "${input}" \
    -o "${output}" \
    -d "${dbsnp_dir}" \
    -ib hg19 \
    -vf "${VARIANT_FILE}" \
    -cf "${chain}" \
    -lt liftOver

# awk '{{if($6<0.5) {{print "chr", $2, "_", $3, "_", $4, "_", $5, "\t", $6, "\t", $10, "\t", $7, "\t", $8, "\t", $9}} else {{print "chr", $2, "_", $3, "_", $4, "_", $5, "\t", 1-$6, "\t", $10, "\t", $7, "\t", $8, "\t", $9}} }}'  OFS="" "${output}_temp" | tail -n +2 | cat <(echo -e "SNP\tMAF\tN\tBETA\tSE\tPVALUE") - > "${output}_temp2"

#     cat <(echo -e "SNP\tMAF\tN\tBETA\tSE\tPVALUE") <(join -1 4 -2 1 <(sort $VARIANT_FILE -k4,4) <(sort $output_temp2 -k1,1) | \
#         awk '{{print $1, $5, $6, $7, $8, $9}}' OFS="\t") > {output}
# OLD METHOD (commented out - used DuckDB + separate liftover):
# # Extract base filename for temporary files
# base_name=$(basename "$input" .tsv.bgz)
# temp_file="${base_name}_temp.tsv"
# 
# echo "Converting GWAS data and performing liftover to GRCh38..."
# 
# # Step 1: Convert using DuckDB (assuming input is GRCh37)
# echo "Step 1: Converting data format..."
# pixi run "duckdb -c \"COPY (SELECT concat('chr', replace(variant, ':', '_')) AS SNP, minor_AF AS MAF, minor_allele AS MINOR_ALLELE, beta AS BETA, se AS SE, n_complete_samples AS N, pval AS PVALUE FROM read_csv('${input}', compression = 'gzip') WHERE NOT low_confidence_variant ) TO '${temp_file}' WITH (FORMAT CSV, DELIMITER E'\\t', HEADER TRUE)\" "
# 
# # Step 2: Perform liftover using R
# echo "Step 2: Performing liftover to GRCh38..."
# pixi run "Rscript liftover.R ${temp_file} ${output}"
# 
# # # Step 3: Clean up temporary file
# echo "Step 3: Cleaning up temporary files..."
# rm -f "${temp_file}"
# 
# echo "Conversion and liftover completed. Output saved to: ${output}"

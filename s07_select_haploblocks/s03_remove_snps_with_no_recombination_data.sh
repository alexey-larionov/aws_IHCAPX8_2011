#!/bin/bash
# s03_remove_snps_with_no_recombination_data.sh

# Alexey Larionov 12Apr2021

# Use:
# ./s03_remove_snps_with_no_recombination_data.sh &> s03_remove_snps_with_no_recombination_data.log

# Stop at runtime errors
set -e

# Start message
echo $0
date
echo ""
echo "Remove ~60k variants with no recombination data"
echo ""

# Files and folders
base_folder="/home/share"

data_folder="${base_folder}/data/s07_select_haploblocks"
input_dataset="${data_folder}/ihcapx8_autosomal_snps"
output_dataset="${data_folder}/ihcapx8_for_ibis"

snps_to_remove="${data_folder}/snps_to_remove.txt"

scripts_folder="${base_folder}/scripts/s07_select_haploblocks"
cd "${scripts_folder}"

# Plink
plink2="${base_folder}/tools/plink/2.0/plink_2.0-a2.3/plink2"

# Progress report
echo "plink2: ${plink2}"
echo "input_dataset: ${input_dataset}"
echo "output_dataset: ${output_dataset}"
echo ""

# Capture ID-s of variants w/o recombination data
echo "Getting ID-s of variants w/o recombination data ..."
awk '$3==0 {print $2}' "${input_dataset}.recomb.bim" > "${snps_to_remove}"
wc -l "${snps_to_remove}"
echo ""

# Filter the plink dataset
echo "Removing variants ..."
"${plink2}" \
  --bfile "${input_dataset}" \
  --exclude "${snps_to_remove}" \
  --make-bed \
  --silent \
  --out "${output_dataset}"
echo ""

# Completion message
echo "Done"
date

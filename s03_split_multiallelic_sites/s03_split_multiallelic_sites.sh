#!/bin/bash
# s01_split_multiallelic_sites.sh

# Alexey Larionov 28Mar2021

# Use:
# ./s01_split_multiallelic_sites.sh &> s01_split_multiallelic_sites.log

# Stop at runtime errors
set -e

# Start message
echo $0
date
echo ""

# Files and folders
base_folder="/home/share"
base_name="IHCAPX8_dragen_joint.hard-filtered.PF.MA-flag"

scripts_folder="${base_folder}/scripts/s03_split_multiallelic_sites"
cd "${scripts_folder}"

data_folder="${base_folder}/data/s03_split_multiallelic_sites"
source_vcf="${data_folder}/${base_name}.vcf.gz"
output_vcf="${data_folder}/${base_name}.MA-split.vcf.gz"

# Progress report
bcftools --version
echo ""
echo "source_vcf: ${source_vcf}"
echo "output_vcf: ${output_vcf}"
echo ""

# Split multiallelic sites
echo "Splitting multiallelic sites ..."
bcftools norm "${source_vcf}" \
-m-any \
-Oz \
-o "${output_vcf}"
echo ""

# Indexing
echo "Indexing ..."
echo ""
bcftools index "${output_vcf}"

# Completion mesage
echo "Done"
date

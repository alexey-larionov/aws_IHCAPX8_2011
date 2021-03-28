#!/bin/bash

# s02_keep_only_pf_variants.sh
# Alexey Larionov, 28Mar2021

# Intended use:
# ./s02_keep_only_pf_variants.sh &> s02_keep_only_pf_variants.log

# Stop at runtime errors
set -e

# Start message
echo $0
date
echo ""

# Files and folders
base_folder="/home/share"
base_name="IHCAPX8_dragen_joint.hard-filtered"

scripts_folder="${base_folder}/scripts/s02_keep_only_pf_variants"
cd "${scripts_folder}"

source_folder="${base_folder}/data/s01_copy_vcf"
source_vcf="${source_folder}/${base_name}.vcf.gz"

output_folder="${base_folder}/data/s02_keep_only_pf_variants"
rm -fr "${output_folder}"
mkdir -p "${output_folder}"
output_vcf="${output_folder}/${base_name}.PF.vcf.gz"

# Progress report
bcftools --version
echo ""
echo "source_vcf: ${source_vcf}"
echo "output_vcf: ${output_vcf}"
echo ""

# Filter vcf
echo "Excluding variants that don't pass all filters ..."
echo ""
bcftools view -Oz -f PASS "${source_vcf}" > "${output_vcf}"

# Index
echo "Indexing ..."
echo ""
bcftools index "${output_vcf}"

# Completion message
echo "Done"
date

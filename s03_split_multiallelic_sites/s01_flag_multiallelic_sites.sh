#!/bin/bash
# s01_flag_multiallelic_sites.sh

# Alexey Larionov 28Mar2021

# Use:
# ./s01_flag_multiallelic_sites.sh &> s01_flag_multiallelic_sites.log

# References & examples

# Stop at runtime errors
set -e

# Start message
echo $0
date
echo ""

# Files and folders
base_folder="/home/share"
base_name="IHCAPX8_dragen_joint.hard-filtered.PF"

scripts_folder="${base_folder}/scripts/s03_split_multiallelic_sites"
cd "${scripts_folder}"

source_folder="${base_folder}/data/s02_keep_only_pf_variants"
source_vcf="${source_folder}/${base_name}.vcf.gz"

output_folder="${base_folder}/data/s03_split_multiallelic_sites"
rm -fr "${output_folder}"
mkdir -p "${output_folder}"
multiallelic_vcf="${output_folder}/${base_name}.MA.vcf.gz"
output_vcf="${output_folder}/${base_name}.MA-flag.vcf.gz"

# Progress report
bcftools --version
echo ""
echo "source_vcf: ${source_vcf}"
echo "multiallelic_vcf: ${multiallelic_vcf}"
echo "output_vcf: ${output_vcf}"
echo ""

# Select multiallelic sites
echo "Selecting multiallelic sites ..."
bcftools view "${source_vcf}" \
--min-alleles 3 \
-Oz \
-o "${multiallelic_vcf}"
echo ""

# Indexing
echo "Indexing ..."
echo ""
bcftools index "${multiallelic_vcf}"

# Flag multiallelic sites
echo "Adding multiallelic flag ..."
bcftools annotate "${source_vcf}" \
--annotations "${multiallelic_vcf}" \
--mark-sites MULTIALLELIC \
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

#!/bin/bash

# s02_explore_vcf.sh
# Alexey Larionov, 28Mar2021

# Intended use:
# ./s02_explore_vcf.sh &> s02_explore_vcf.log

# Stop at runtime errors
set -e

# Start message
echo $0
date
echo ""

# Files and folders
base_folder="/home/share"
base_name="IHCAPX8_dragen_joint.hard-filtered.PF"

scripts_folder="${base_folder}/scripts/s02_keep_only_pf_variants"
cd "${scripts_folder}"

data_folder="${base_folder}/data/s02_keep_only_pf_variants"

vcf="${data_folder}/${base_name}.vcf.gz"

stats_folder="${data_folder}/bcfstats"
mkdir -p "${stats_folder}"
stats_file="${stats_folder}/${base_name}.vchk"

annotations="${data_folder}/${base_name}.annotations.txt"
summary="${data_folder}/${base_name}.summary.txt"

# Progress report
bcftools --version
echo ""
echo "vcf: ${vcf}"
echo "data_folder: ${data_folder}"
echo ""

# Count variants
echo "Variant counts:"
echo ""
bcftools +counts "${vcf}"
echo ""

# Calculate bcfstats
echo "Calculating bcfstats..."
echo ""
bcftools stats -s- "${vcf}" > "${stats_file}"

# Plot bcfstats
# requires python3 with matplotlib (for plotting) and texlive (for making PDF)
echo "Making bcfstats plots..."
echo ""
plot-vcfstats -p "${stats_folder}" "${stats_file}"
echo ""

# Extract annotations
echo "Extracting annotations ..."
echo ""
bcftools query -H -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER\t%QUAL\t%FORMAT\n' "${vcf}" > "${annotations}"

echo "--- FILTER ---"
echo ""
awk 'NR>1 {print $6}' "${annotations}" | sort |  uniq -c | sort -nr
echo ""

# Make summary with vcf header and selected records
echo "Making summary ..."
echo ""
bcftools view -h "${vcf}" > "${summary}"
bcftools view -H "${vcf}" | head >> "${summary}"
echo "..." >> "${summary}"
bcftools view -H "${vcf}" | tail >> "${summary}"

# Completion message
echo "Done"
date
echo ""

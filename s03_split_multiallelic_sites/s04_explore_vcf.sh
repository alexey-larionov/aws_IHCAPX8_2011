#!/bin/bash

# s04_explore_vcf.sh
# Alexey Larionov, 28Mar2021

# Intended use:
# ./s04_explore_vcf.sh &> s04_explore_vcf.log

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
split_vcf="${data_folder}/${base_name}.MA-split.vcf.gz"

stats_folder="${data_folder}/bcfstats"
rm -fr "${stats_folder}"
mkdir -p "${stats_folder}"
stats_file="${stats_folder}/${base_name}.MA-split.vchk"

annotations="${data_folder}/${base_name}.MA-split.annotations.txt"
summary="${data_folder}/${base_name}.MA-split.summary.txt"

# Progress report
bcftools --version
echo ""
echo "source_vcf: ${source_vcf}"
echo "split_vcf: ${split_vcf}"
echo "data_folder: ${data_folder}"
echo ""

# Count variants

echo "Variant counts in source VCF:"
echo ""
bcftools +counts "${source_vcf}"
echo ""

echo "Variant counts in split VCF:"
echo ""
bcftools +counts "${split_vcf}"
echo ""

# Calculate bcfstats in split VCF
echo "Calculating bcfstats in split VCF..."
bcftools stats -s- "${split_vcf}" > "${stats_file}"
echo ""

# Plot bcfstats
# requires python3 with matplotlib (for plotting) and texlive (for making PDF)
echo "Making bcfstats plots..."
echo ""
plot-vcfstats -p "${stats_folder}" "${stats_file}"
echo ""

# Extract annotations

echo "Extracting annotations from split VCF ..."
bcftools query -H -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER\t%QUAL\t%FORMAT\n' "${split_vcf}" > "${annotations}"
echo ""

# Look at specific examples

echo "A bi-allelic site in split vcf"
echo ""
bcftools view -H -r chr12:128667912 "${split_vcf}"
echo ""

echo "A multiallelic site in source vcf"
echo ""
bcftools view -H -r chr1:15274 "${source_vcf}"
echo ""

echo "A multiallelic site in split vcf"
echo ""
bcftools view -H -r chr1:15274 "${split_vcf}"
echo ""

echo "A warning multiallelic site in source vcf"
echo ""
bcftools view -H -r chr1:12866791 "${source_vcf}"
echo ""

echo "A warning multiallelic site in split vcf"
echo ""
bcftools view -H -r chr1:12866791 "${split_vcf}"
echo ""

# Make summary with vcf header and selected records in split VCF
echo "Making summary of split VCF ..."
bcftools view -h "${split_vcf}" > "${summary}"
bcftools view -H "${split_vcf}" | head >> "${summary}"
echo "..." >> "${summary}"
bcftools view -H "${split_vcf}" | tail >> "${summary}"
echo ""

# Completion message
echo "Done"
date
echo ""

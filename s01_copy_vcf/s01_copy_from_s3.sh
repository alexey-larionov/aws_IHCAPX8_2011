#!/bin/bash

# s01_copy_from_s3.sh
# Alexey Larionov, 25Mar2021

# Intended use:
# ./s01_copy_from_s3.sh &> s01_copy_from_s3.log

# Data transfer speed ~60MiB/sec, total transfer time for 1GB VCF < 1 min

# Stop at runtime errors
set -e

# Start message
echo "Copy VCF from s3"
date
echo ""

# Folders
base_folder="/home/share"
data_folder="${base_folder}/data/s01_vcf"
mkdir -p "${data_folder}"

# Copy data
cd "${data_folder}"
aws s3 cp "s3://ihcapx8/vcfs/IHCAPX8_dragen_joint.hard-filtered.vcf.gz" "${data_folder}/" --quiet
aws s3 cp "s3://ihcapx8/vcfs/IHCAPX8_dragen_joint.hard-filtered.vcf.gz.md5sum" "${data_folder}/" --quiet
aws s3 cp "s3://ihcapx8/vcfs/IHCAPX8_dragen_joint.hard-filtered.vcf.gz.tbi" "${data_folder}/" --quiet

# Check vcf file
md5sum IHCAPX8_dragen_joint.hard-filtered.vcf.gz
cat IHCAPX8_dragen_joint.hard-filtered.vcf.gz.md5sum
echo ""

# Completion message
echo "Done"
date

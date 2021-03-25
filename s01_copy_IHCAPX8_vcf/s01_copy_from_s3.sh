#!/bin/bash

# s01_copy_from_s3.sh
# Alexey Larionov, 24Mar2021

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
scripts_folder="${base_folder}/scripts/s01_copy_data"
cd "${scripts_folder}"

# Copy data
aws s3 cp "s3://ihcapx8/vcfs/IHCAPX8_dragen_joint.hard-filtered.vcf.gz" "${data_folder}/"
aws s3 cp "s3://ihcapx8/vcfs/IHCAPX8_dragen_joint.hard-filtered.vcf.gz.md5sum" "${data_folder}/"
aws s3 cp "s3://ihcapx8/vcfs/IHCAPX8_dragen_joint.hard-filtered.vcf.gz.tbi" "${data_folder}/"

# Completion message
echo "Done"
date

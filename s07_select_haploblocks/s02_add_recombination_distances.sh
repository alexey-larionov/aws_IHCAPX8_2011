#!/bin/bash
# s02_add_recombination_distances.sh

# Alexey Larionov 12Apr2021

# Use:
# ./s02_add_recombination_distances.sh &> s02_add_recombination_distances.log

# Stop at runtime errors
set -e

# Start message
echo $0
date
echo ""
echo "Adding recombination distances"
echo ""

# Files and folders
base_folder="/home/share"

data_folder="${base_folder}/data/s07_select_haploblocks"
plink_dataset="${data_folder}/ihcapx8_autosomal_snps"

scripts_folder="${base_folder}/scripts/s07_select_haploblocks"
cd "${scripts_folder}"

ibis_folder="${base_folder}/tools/ibis"
map_folder="/home/share/resources/recombination_maps/b38/recomb-hg38"

# Progress report
echo "plink_dataset: ${plink_dataset}"
echo "ibis_folder: ${ibis_folder}"
echo "map_folder: ${map_folder}"
echo ""

# Update bim
echo "Updating bim ..."
"${ibis_folder}/add-map-plink.pl" \
  "${plink_dataset}.bim" \
  "${map_folder}/genetic_map_GRCh38_merged_no_chrX.tab" \
  > "${plink_dataset}.recomb.bim"
echo ""

# Completion message
echo "Done"
date

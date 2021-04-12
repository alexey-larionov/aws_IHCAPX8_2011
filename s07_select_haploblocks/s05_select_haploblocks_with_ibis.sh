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
echo "Import VCF to PLINK (bed-bim-fam file-set)"
echo ""

# Files and folders
base_folder="/home/share"

data_folder="${base_folder}/data/s07_select_haploblocks"
plink_dataset="${data_folder}/ihcapx8"

scripts_folder="${base_folder}/scripts/s07_select_haploblocks"
cd "${scripts_folder}"

map_directory="/home/share/resources/recombination_maps/GRCh38"

# IBIS
ibis="${base_folder}/tools/ibis/ibis"


# Progress report
echo "source_vcf: ${source_vcf}"
echo "plink_dataset: ${plink_dataset}"
echo "plink: ${plink2}"
echo ""

# Update bim
./add-map-plink.pl "${plink_dataset}.bim" ${map_directory}/plink.chr{1..22}.GRCh38.map > "${plink_dataset}.new.bim"
genetic_map_GRCh37_chr{1..22}.txt > new.bim

# Completion message
echo "Done"
date

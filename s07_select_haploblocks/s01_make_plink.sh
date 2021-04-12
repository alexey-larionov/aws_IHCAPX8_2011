#!/bin/bash
# s01_make_plink.sh

# Alexey Larionov 12Apr2021

# Use:
# ./s01_make_plink.sh &> s01_make_plink.log

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

source_folder="${base_folder}/data/s04_annotate"
source_vcf="${source_folder}/IHCAPX8_dragen_joint.hard-filtered.PF.MA-flag.MA-split.ID.ClinVar.std-Chr.Reheaded.VEP.split-VEP.vcf.gz"

plink_dataset="${base_folder}/data/s07_select_haploblocks/ihcapx8_autosomal_snps"

scripts_folder="${base_folder}/scripts/s07_select_haploblocks"
cd "${scripts_folder}"

# Plink
plink2="${base_folder}/tools/plink/2.0/plink_2.0-a2.3/plink2"

# Progress report
echo "source_vcf: ${source_vcf}"
echo "plink_dataset: ${plink_dataset}"
echo "plink: ${plink2}"
echo ""

# Import VCF to PLINK
# --vcf-half-call describes what to do with genotypes like 0/.
# --allow-no-sex suppresses warning about missed sex
# --double-id puts sample name to both Family-ID and Participant-ID
# --autosome excludes all unplaced and non-autosomal variants
# --silent suppresses very verbous ouput to the "out" file (log file is still avaialble in the data folder)
echo "Making plink dataset ..."
"${plink2}" \
  --vcf "${source_vcf}" \
  --vcf-half-call "missing" \
  --double-id \
  --allow-no-sex \
  --autosome \
  --snps-only \
  --make-bed \
  --silent \
  --out "${plink_dataset}"

# Completion message
echo ""
echo "Done"
date

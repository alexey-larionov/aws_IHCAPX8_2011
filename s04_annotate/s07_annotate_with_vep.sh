#!/bin/bash
# s07_annotate_with_vep.sh

# Alexey Larionov 10Apr2021

# Intended use:
# ./s07_annotate_with_vep.sh &> s07_annotate_with_vep.log

# Notes

# - Use 4xLarge instance with 16 cores to run within 10-12 hrs
# - Dont edit script during its execution (somehow it may confuse bash)

# For tests
# - Use an instance with at least 4 cores and 16GB RAM (i.e use xlarge and above)
# - The VEP annotation takes 40-50 min per VCF with ~250k variants using 4 threads
#   (on m4.xlarge instance)

# Mount volume with CADD data before running:
# 1) Make EBS volume from snapshot
#    - gp2 (default) type, in the same availability zone as instance,
#    - optionally add tag Name=CADD (just for convenience)
# 2) Attach volume to the running instance
#    - using AWS console for EC2-Volume
# 3) Mount the volume using instance terminal
#    - lsblk # Check that the voume is attached
#    - ls /home/share/resources/CADD # Check that the mount point is available
#    - sudo mount /dev/xvdf /home/share/resources/CADD # Mount the volume
#    - ls /home/share/resources/CADD -all # Check that volume has been mounted (and see the folder permissions)
# 4) Use screen to run the script (because it takes hours to run)
#    - screen -ls
#    - screen -S vep
#    - ./s07_annotate_with_vep.sh &> s07_annotate_with_vep.log
#    - ctr-a d
#    - screen -r vep
#    - exit

# Unmount and terminate the volume after running:
# 1) In the instance terminal: sudo umount /home/share/resources/CADD
# 2) In the AWS console for EC2-Volume: Detach, then Delete
# 3) Do NOT make AMI (Image) before detaching CADD volume (otherwise CADD volume will be included in the AMI & snapshot)

# Stop at runtime errors
set -e

# Start message
echo $0
date
echo ""

# Files and folders
base_folder="/home/share"
base_name="IHCAPX8_dragen_joint.hard-filtered.PF.MA-flag.MA-split.ID.ClinVar.std-Chr.Reheaded"

scripts_folder="${base_folder}/scripts/s04_annotate"
cd "${scripts_folder}"

data_folder="${base_folder}/data/s04_annotate"

source_vcf="${data_folder}/${base_name}.vcf.gz"
output_vcf="${data_folder}/${base_name}.VEP.vcf.gz"
vep_report="${data_folder}/${base_name}.VEP.html"

vep_cache_info="${data_folder}/vep_cache_info.txt"

# VEP script
vep="${base_folder}/tools/ensembl/vep-102-b38/vep"

# VEP cache
# cache_folder isn't needed because cache was installed to the default location
# i.e cache_folder="/home/ec2-user/.vep
cache_version="102"
cache_assembly="GRCh38"

# VEP plugins folder
# plugins_folder isn't needed because plugins were installed to the default location
# i.e. plugins_folder="/home/ec2-user/.vep"

# Data for plugins
# Location should correspond to actually mounted data
cadd_data_folder="${base_folder}/resources/CADD/GRCh38/v_1.6"
cadd_snv_data="${cadd_data_folder}/whole_genome_SNVs.tsv.gz"
cadd_indels_data="${cadd_data_folder}/gnomad.genomes.r3.0.indel.tsv.gz"

# Reference genome
# A reference genome is provided in cache, i.e.
# b38_fasta="/home/ec2-user/.vep/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
# Most likely, VEP uses this reference by default (could be clarifyed in docs/ asking)
# The exact reference used by Illumina BaseSpace Dragen is not yet known:
# so it may be useful to ask Illumina what reference genome could be used

# Num of threads (should fit the hardware)
n_threads="14"

# Progress report
echo "--- Files ---"
echo ""
echo "source_vcf: ${source_vcf}"
echo "output_vcf: ${output_vcf}"
echo "vep_report: ${vep_report}"
echo ""
echo "--- VEP cache ---"
echo ""
#echo "cache_folder: ${cache_folder}" # not needed because cache is in the default location in ~/.vep
echo "cache_version: ${cache_version}"
echo "cache_assembly: ${cache_assembly}"
echo ""
"${vep}" --show_cache_info
echo ""
cp "/home/ec2-user/.vep/homo_sapiens/102_GRCh38/info.txt" "${vep_cache_info}"
echo "See more information about used cache in the following file:"
echo "${vep_cache_info}"
echo ""
echo "--- VEP plugins ---"
echo ""
#echo "plugins_folder: ${plugins_folder}" # not needed because cache is in the default location in ~/.vep
#echo ""
echo "CADD annotation files:"
echo "${cadd_snv_data}"
echo "${cadd_indels_data}"
echo ""
#echo "--- Reference ---"
#echo ""
#echo "b38_fasta: ${b38_fasta}" # not needed because the default cache reference is used
#echo ""
echo "--- Other settings ---"
echo ""
echo "n_threads: ${n_threads}"
echo ""

# Annotate VCF
echo "Annotating ..."
"${vep}" \
--input_file "${source_vcf}" \
--output_file "${output_vcf}" \
--stats_file "${vep_report}" \
--vcf \
--force_overwrite \
--compress_output bgzip \
--fork "${n_threads}" \
--offline \
--cache \
--species homo_sapiens \
--cache_version "${cache_version}" \
--assembly "${cache_assembly}" \
--pick \
--gencode_basic \
--check_ref \
--everything \
--nearest symbol \
--total_length \
--check_existing \
--exclude_null_alleles \
--plugin CADD,"${cadd_snv_data}","${cadd_indels_data}"

echo ""

# Not needed because the default locations are used:
#  --dir_cache "${cache_folder}" \
#  --dir_plugins "${plugins_folder}" \

# Index annotated vcf
echo "Indexing ..."
bcftools index "${output_vcf}"
echo ""

# Added VEP annotations
echo "Added VEP annotations:"
echo ""
bcftools +split-vep -l "${output_vcf}"
echo ""

# Completion message
echo "Done"
date

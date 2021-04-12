#!/bin/bash
# s01_render_rmd.sh

# Alexey Larionov 10Apr2021

# Intended use:
# This script shoudl be run in screen as r-user 
# screen -S rnd
# su - r-user # aws-2021
# ./s01_render_rmd.sh &> s01_render_rmd.log

# Reference:
# https://stackoverflow.com/questions/35077224/paths-in-ssh-versus-rstudio-server-for-pandoc-and-knitr

# Stop at runtime errors
set -e

# Start message
echo $0
date
echo ""

# Render rmd
PATH=/usr/lib/rstudio-server/bin/pandoc/:$PATH
cd "/home/share/scripts/s05_read_to_R"
Rscript -e "rmarkdown::render('s01_import_vcf_to_R.Rmd')"

# Completion message
echo ""
echo "Done"
date

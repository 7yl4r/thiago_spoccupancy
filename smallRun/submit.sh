#!/bin/bash
module load apps/R/4.3.0
module load apps/singularity/3.5.0

#apt-get install -y cmake libssl-dev libudunits2-dev libgdal-dev # might also need `gdal-bin`
yum install -y cmake openssl-devel udunits2-devel gdal-devel

# Install R packages inside the Singularity container
singularity exec /apps/R/4.3.0/R_4.3.0.sif Rscript -e "
  install.packages(
    c('spOccupancy', 'stars', 'ggplot2'),
    repos='https://cloud.r-project.org/'
  )
"

#Rscript -e "install.packages('spOccupancy')"
#Rscript -e "install.packages('stars')"
#Rscript -e "install.packages('ggplot2')"

singularity exec /apps/R/4.3.0/R_4.3.0.sif Rscript 5.Analysis_spOccupancy_MultiSp_SpatInteg_Summer.R

## setup
starting from ubuntu 24 with r-base installed
* `sudo apt install -y cmake libssl-dev libudunits2-dev libgdal-dev`
   * might also need `gdal-bin`
* `Rscript -e "install.packages("spOccupancy")`
* `Rscript -e "install.packages("stars")`
* `Rscript -e "install.packages("ggplot2")`

## run
* `Rscript 5...`5.Analysis_spOccupancy_MultiSp_SpatInteg_Summer.R

### SLURM

#### set up new run
`mkdir` run##, `cd` into it, then:

```bash
mkdir old-slurm
touch README.md

ln -s ../data/Shapes Shapes
ln -s ../data/spOccupancy_MultiSpp_FullArea spOccupancy_MultiSpp_FullArea
ln -s ../data/Grid_OccEnv_Seasonal.txt Grid_OccEnv_Seasonal.txt

cp ../run03/5.Analysis_spOccupancy_MultiSp_SpatInteg_Winter.R .
cp ../run03/5.Analysis_spOccupancy_MultiSp_SpatInteg_Summer.R .
cp ../run03/submit.sh .
```
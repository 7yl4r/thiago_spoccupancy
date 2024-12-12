SDM modeling built from acoustic tag data of 69 species in the Florida region.
R code packaged using Singularity to run on USF's SLURM supercomputer CIRCE.
This is a collaborative project under the SE US MBON. 

This repo includes source code under `./src/` and details on each supercomputer run taken in the debugging process (`./r*/` dirs).   

## setup
starting from ubuntu 24 with r-base installed
* `sudo apt install -y cmake libssl-dev libudunits2-dev libgdal-dev`
   * might also need `gdal-bin`
* `Rscript -e "install.packages("spOccupancy")`
* `Rscript -e "install.packages("stars")`
* `Rscript -e "install.packages("ggplot2")`

## run
* `Rscript 5...`5.Analysis_spOccupancy_MultiSp_SpatInteg_Summer.R

### Singularity
Singularity is used to package the application for a SLURM supercomputer.

```bash
# Build `run.sif` from `Singularity.def` on local machine with sudo:
sudo singularity build run.sif Singularity.def

# transfer .sif to supercomputer
rsync -hazv run.sif tylarmurray@circe.rc.usf.edu:.


```

### SLURM

#### set up new run
`mkdir` run##, `cd` into it, then:

```bash
TEMPLATE=run03
mkdir old-slurm
cat "new run from template $TEMPLATE" > README.md

ln -s ../data/Shapes Shapes
ln -s ../data/spOccupancy_MultiSpp_FullArea spOccupancy_MultiSpp_FullArea
ln -s ../data/Grid_OccEnv_Seasonal.txt Grid_OccEnv_Seasonal.txt

cp ../src/5.Analysis_spOccupancy_MultiSp_SpatInteg.R .

cp ../$TEMPLATE/submit.sh .

git add README.md 5.Analysis_spOccupancy_MultiSp_SpatInteg.R submit.sh
git add -f Grid_OccEnv_Seasonal.txt Shapes spOccupancy_MultiSpp_FullArea
```

NOTE: for *winter runs, the file is still named `*_Summer.R`, just use a different `TEMPLATE`.

# Authorship
* Thiago Belisario D'Araujo Couto
* Viviane Zulian
* Tylar Murray
* Neil Hammerschlag

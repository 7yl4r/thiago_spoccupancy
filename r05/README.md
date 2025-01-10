This run has some issue with the revised R code resulting in the following error:

```
 [ reached getOption("max.print") -- omitted 15342 rows ]
[1] "sf"         "data.frame"
Deleting layer `SDM_Shape_winter_NA' using driver `ESRI Shapefile'
Writing layer `SDM_Shape_winter_NA' to data source 
  `SDM_Shape_winter_NA.shp' using driver `ESRI Shapefile'
Writing 27841 features with 3 fields and geometry type Multi Polygon.
Warning message:
In abbreviate_shapefile_names(obj) :
  Field names abbreviated for ESRI Shapefile driver
[1] "1"
```

The summer+winter merged script is having some issue where all species names are being set to NA.
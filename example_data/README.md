# SkyStrat San Rafael Swell demo dataset 

### Alex Tye, 5/29/2026

### Contents

The files are housed within two .zip files:

**gdb.zip** - zip file of Skystrat_packagetoshare.gdb, containing script inputs and outputs
- **SRswell_studyarea** - polygon feature class of the study area, a required INPUT
- **strike_dip** - point feature class with strike and dip measurements across the study area, a required INPUT
- **strike_dip_withheight** - point feature class with strike/dip data and inferred stratigraphic heights, an optional OUTPUT

**rasters.zip** - zip file containing input and output rasters for the script
- **output_USGS1m.tif** - USGS 1m DEM of the study area (United States Geological Survey, 2021), a required INPUT
- **wedge.tif** - a GeoTiff mapping which circle sector each point on the DEM falls within, an optional OUTPUT
- **height.tif** - a GeoTiff of stratigraphic height, an optional OUTPUT
- **testoutput.pdf** - not a raster, but a PDF showing the profile view of folding, a mandatory OUTPUT


### References

United States Geological Survey (2021). United States Geological Survey 3D Elevation Program 1 meter Digital Elevation Model. Distributed by OpenTopography. https://doi.org/10.5069/G9NP22NT. Accessed 2026-05-29

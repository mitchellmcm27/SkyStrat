# SkyStrat San Rafael Swell demo dataset 

### Alex Tye, 5/29/2026

### Description

The files in this folder are a demonstration dataset that show the use of the SkyStrat tool, which uses the Busk (1929) cross section construction method to calculate and map stratigraphic height across a landscape. The files include a GeoDataBase (.gdb) file inside a zip folder, which contains three vector datasets for the example problem, a polygon of the study area and a set of strike and dip measurements, which are both required inputs for the tool. Also contained is a copy of the strike and dip measurement vector dataset that has stratigraphic height calculated and listed for each strike and dip measurement station, which is an optional output from the script.

In addition to the vector datasets, a second zip folder contains three raster files (GeoTiff) for the example problem. A digital elevation model (DEM) is included, which is a required input. Also shown are two output rasters, one of which shows the identity of the wedge or circle sector that each DEM point was assigned to in the analysis, and the second is a raster of stratigraphic height calculated using the tool. Finally, a PDF is included which is a figure created and output by the tool.

If the user runs the `BuskMethod` tool using the files marked INPUT, an output that is identical to those marked OUTPUT should be created.

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

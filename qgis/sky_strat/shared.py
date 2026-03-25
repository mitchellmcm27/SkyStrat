
def _sample_raster_at_xy(ds, x, y):
    """
    Sample a single cell value from an open GDAL dataset at geographic (x, y).
    Returns None if the coordinate is outside the raster extent or on a NoData cell.
    Replaces arcpy.management.GetCellValue / arcpy.RasterToNumPyArray per-cell lookup.
    """
    gt = ds.GetGeoTransform()
    
    col = int((x - gt[0]) / gt[1])
    row = int((y - gt[3]) / gt[5])       # gt[5] < 0 for north-up rasters
    ncols = ds.RasterXSize
    nrows = ds.RasterYSize
    if col < 0 or col >= ncols or row < 0 or row >= nrows:
        return None
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray(col, row, 1, 1)
    if arr is None:
        return None
    v = float(arr[0, 0])
    nodata = band.GetNoDataValue()
    if nodata is not None and v == nodata:
        return None
    return v
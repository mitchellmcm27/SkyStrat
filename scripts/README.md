# BuskMethod.pyt — Busk Method Toolbox

An ArcGIS Python Toolbox (`.pyt`) that implements the **Busk Method** for calculating stratigraphic heights from strike/dip measurements using cylindrical fold geometry. Given bedding orientation measurements distributed across a study area, this tool determines the stratigraphic height (vertical position in the rock sequence) of every point and can produce stratigraphic height rasters for the entire study area.

---

## Table of Contents

1. [Overview](#overview)
2. [Dependencies](#dependencies)
3. [Tool: Busk Downplunge View](#tool-busk-downplunge-view)
   - [Input Parameters](#input-parameters)
   - [Output Parameters](#output-parameters)
4. [Workflow](#workflow)
5. [Key Algorithms](#key-algorithms)
   - [Fold Axis Calculation](#1-fold-axis-calculation-eigenvalue-method)
   - [Downplunge Projection](#2-downplunge-projection-onto-profile-plane)
   - [Wedge Analysis](#3-wedge-analysis)
   - [Stratigraphic Height Propagation](#4-stratigraphic-height-propagation)
   - [DEM Cell Assignment](#5-dem-cell-assignment)
6. [Outputs](#outputs)
7. [Data Constraints and Assumptions](#data-constraints-and-assumptions)
8. [Helper Functions Reference](#helper-functions-reference)

---

## Overview

The Busk Method is a classical structural geology technique for constructing cross-sections through folded strata using a concentric (cylindrical) fold model. This toolbox automates that process within the ArcGIS environment:

- Projects field strike/dip measurements onto a **downplunge profile plane** perpendicular to the fold axis
- Uses the resulting 2D geometry to propagate **stratigraphic heights** through adjacent measurement pairs (wedges)
- Assigns stratigraphic heights to every **DEM cell** within the study area
- Outputs a **PDF visualization** of the downplunge projection along with optional rasters and feature classes

---

## Dependencies

| Library | Purpose |
|---|---|
| `arcpy` | ArcGIS spatial operations, raster/feature I/O |
| `numpy` | Array math, vector operations, eigenvalue decomposition |
| `matplotlib` | Downplunge view visualization |
| `matplotlib.backends.backend_pdf.PdfPages` | PDF export |
| `os`, `tempfile` | File system operations |

Requires an **ArcGIS Pro** installation with a Spatial Analyst license (for `ExtractByMask` and `GetCellValue`).

---

## Tool: Busk Downplunge View

**Label**: `Busk Downplunge View`
**Class**: `BuskDownplungeView`
**Description**: Creates a downplunge projection of strike/dip measurements and optionally computes stratigraphic heights for measurement points and DEM cells.

### Input Parameters

| # | Name | Type | Required | Description |
|---|---|---|---|---|
| 0 | Strike and Dip Point Feature Class | Point Feature Layer | Yes | Input point features containing strike/dip measurements |
| 1 | Strike Field | Field | Yes | Field with strike values in azimuth degrees (0–360°) |
| 2 | Dip Field | Field | Yes | Field with dip values in degrees (0–90°) |
| 3 | Overturned Field | Field | No | Field flagging overturned beds (`1` = overturned); flips the stratigraphic height vector |
| 4 | Label Field | Field | Yes | Field used to identify each measurement point (defaults to OBJECTID) |
| 5 | Stratigraphic Height Field | Field | No | Field with a known stratigraphic height for **one** reference measurement; used to calibrate all calculated heights |
| 6 | DEM Raster | Raster Layer | Yes | Digital Elevation Model for elevation extraction and raster output |
| 7 | Study Area | Polygon Feature Layer | Yes | Boundary polygon clipping the analysis region |
| 8 | Fold Trend | Double (0–360°) | No | Fold axis azimuth; auto-calculated from poles if not provided |
| 9 | Fold Plunge | Double (0–90°) | No | Fold axis plunge angle; auto-calculated from poles if not provided |
| 10 | Plot DEM Cells | Boolean | No | If `True`, overlays color-coded DEM cell centroids on the profile plot (default: `False`) |

### Output Parameters

| # | Name | Type | Required | Description |
|---|---|---|---|---|
| 10 | Output PDF | File (.pdf) | Yes | Downplunge projection visualization |
| 11 | Wedge Raster Output | Raster Dataset | No | Integer raster showing wedge/rectangle index assignment per DEM cell |
| 12 | Output Strike/Dip Feature Class | Feature Class | No | Copy of input features with appended `strat_height_out` field |
| 13 | Stratigraphic Height Raster Output | Raster Dataset | No | Float32 raster of calculated stratigraphic heights for all DEM cells |

---

## Workflow

```
 1. Validate inputs (check for projected/UTM coordinate systems)
 2. Read strike/dip data; extract missing elevations from DEM
 3. Calculate fold axis (from poles eigenvalue method) or accept user-provided values
 4. Project all measurement points onto the downplunge profile plane
 5. Analyze wedges between adjacent measurement pairs
 6. [If raster output requested] Extract DEM cell centroids within study area
 7. Find intersection points for non-parallel wedge pairs
 8. Propagate stratigraphic heights across measurement points (left → right in profile)
 9. [If raster output requested] Assign each DEM cell to its wedge/rectangle
10. [If raster output requested] Calculate stratigraphic height for each DEM cell
11. [If requested] Write wedge assignment raster
12. [If requested] Write stratigraphic height raster
13. [If requested] Write output feature class with calculated heights
14. Render and save PDF downplunge view visualization
15. Report final statistics via arcpy messages
```

---

## Key Algorithms

### 1. Fold Axis Calculation (Eigenvalue Method)

Assumes **cylindrical fold geometry**: bedding poles form a great circle in orientation space, and the fold axis is the normal to that great circle.

- Converts each strike/dip measurement to its pole unit vector (East, North, Up)
- Constructs the covariance matrix of all poles
- Performs eigenvalue decomposition; the **eigenvector with the smallest eigenvalue** is the fold axis
- Converts the result back to trend/plunge format

### 2. Downplunge Projection onto Profile Plane

Defines an orthonormal 2D coordinate system on the profile plane (perpendicular to the fold axis):

- `profile_x`: horizontal axis in the profile plane
- `profile_y`: vertical axis in the profile plane

Each 3D measurement point is projected via dot products with these basis vectors. For each measurement, the bedding plane is intersected with the profile plane to obtain the **bedding trace direction** and a **stratigraphic height vector** (perpendicular to bedding, pointing up-section for upright beds and down-section for overturned beds).

### 3. Wedge Analysis

Measurement points are ordered by their `profile_x` coordinate (left to right). Adjacent pairs define **wedges** or **rectangles**:

- **Rectangle (Parallel Beds)**: Bedding attitude vectors differ by < 0.01°. Stratigraphic height changes linearly by projecting the displacement vector onto the stratigraphic height vector.
- **Wedge (Intersecting)**: Attitude vectors are not parallel. The intersection of the two stratigraphic height lines is calculated using Cramer's rule. Wedge type is classified by the sign of the line parameters at intersection:
  - **UP-UP**: Both parameters positive — heights decrease with distance from intersection
  - **DOWN-DOWN**: Both parameters negative — heights increase with distance from intersection
  - **Invalid**: Opposite signs — geometrically inconsistent pair, skipped

### 4. Stratigraphic Height Propagation

Requires exactly **one** measurement with a known stratigraphic height (or heights are relative, starting at 0 for the leftmost point).

**Parallel (Rectangle)**:
```
height_right = height_left + dot(displacement_vector, strat_height_vector)
```

**Intersecting Wedge — UP-UP**:
```
height = -distance_to_intersection + constant
constant = height_left + distance(left_point, intersection)
```

**Intersecting Wedge — DOWN-DOWN**:
```
height = +distance_to_intersection + constant
constant = height_left - distance(left_point, intersection)
```

A single **correction constant** is applied at the end to shift all values to match the known reference height.

### 5. DEM Cell Assignment

All DEM grid cell centroids are tested against each wedge/rectangle:

**Rectangle containment test**:
- Translate to left-point origin
- Rotate so the stratigraphic height vector aligns with +Y
- Check: `0 ≤ x ≤ right_x` and `0 ≤ y ≤ right_y`

**Wedge containment test**:
- Translate to intersection point origin
- Check angular containment using dot products against boundary vectors toward each measurement point

Cells falling in multiple regions are flagged as **ambiguous** (value `-2`); unassigned cells receive `-1`.

---

## Outputs

### PDF Visualization

A single-page PDF showing the downplunge projection with:

| Element | Style |
|---|---|
| DEM cell centroids (optional) | Colored scatter points by wedge index |
| Measurement points | Red dots with black borders |
| Bedding traces | Blue lines |
| Stratigraphic height vectors | Gold/yellow lines |
| Rectangle boundaries | Green outlines |
| Wedge boundaries | Black dashed lines to intersection point |
| Point labels | Text from label field |
| Scale bar | Dynamic round-number scale (lower-right) |
| Fold axis info | Shown in legend as `Fold Axis: {trend}° / {plunge}°` |

The X-axis is inverted following the conventional downplunge view orientation.

### Wedge Assignment Raster

- **Type**: Int16
- **NoData**: `-9999`
- **Values**: Wedge/rectangle index (0-based) for each DEM cell
- Matches the DEM spatial reference and cell size, clipped to the study area

### Stratigraphic Height Raster

- **Type**: Float32
- **NoData**: `-9999`
- **Values**: Calculated stratigraphic height in meters
- Same spatial reference, cell size, and extent as the wedge raster

### Output Strike/Dip Feature Class

- Copy of the input feature class
- Adds field `strat_height_out` (Double, alias: "Calculated Stratigraphic Height")
- Populated using the label field as the join key

---

## Data Constraints and Assumptions

1. **Projected coordinates required**: All inputs (feature classes, DEM) must use a projected coordinate system. UTM is recommended and verified; a warning is issued for non-UTM projections.
2. **Units in meters**: All distance and elevation values are assumed to be in meters.
3. **Cylindrical fold geometry**: The tool assumes a single consistent fold axis across the study area.
4. **One reference height**: At most one measurement may carry a known stratigraphic height for calibration.
5. **Complete DEM coverage**: NoData gaps in the DEM within the study area will result in unassigned cells in the output rasters.
6. **Adjacent pair ordering**: Stratigraphic heights are propagated through pairs ordered by profile X position; a connected chain of overlapping wedges is required for full coverage.

---

## Helper Functions Reference

| Function | Description |
|---|---|
| `read_strike_dip_data()` | Reads all measurement attributes and validates inputs |
| `extract_z_from_dem()` | Queries DEM raster for elevation at point locations |
| `calculate_fold_axis()` | Eigenvalue-based fold axis from poles |
| `calculate_projected_attitudes()` | Projects measurements onto profile plane |
| `analyze_wedges()` | Classifies adjacent pairs as parallel or intersecting |
| `extract_dem_cell_centroids()` | Extracts all DEM cell centroids within study area |
| `calculate_stratigraphic_heights_in_profile()` | Propagates heights through measurement chain |
| `assign_dem_cells_to_wedges()` | Assigns each DEM cell to a wedge/rectangle |
| `calculate_dem_stratigraphic_heights()` | Computes height for each assigned DEM cell |
| `create_output_strike_dip_fc()` | Writes output feature class with calculated heights |
| `create_wedge_assignment_raster()` | Writes integer wedge-index raster |
| `create_strat_height_raster()` | Writes float stratigraphic height raster |
| `create_downplunge_view()` | Renders and saves the PDF visualization |
| `find_line_intersection()` | 2D line intersection via Cramer's rule |
| `find_line_intersection_with_params()` | Same, also returns t₁/t₂ parameters for wedge classification |
| `check_utm_projection()` | Validates projected coordinate system on inputs |
| `add_scale_bar()` | Adds a dynamic scale bar to the matplotlib plot |

---

## ArcGIS Dependency Analysis

The table below classifies every method by how tightly it is coupled to `arcpy`. Functions with no arcpy dependency can be extracted and reused in other GIS environments (QGIS/PyQGIS, GDAL/OGR, GeoPandas, standalone scripts, etc.) with little or no modification.

### Fully arcpy-dependent — ArcGIS required

These functions rely on arcpy for core data I/O, spatial operations, or tool-interface registration. They cannot run outside ArcGIS without replacing those arcpy calls with equivalents from another library.

| Function | arcpy APIs used | What would need replacing |
|---|---|---|
| `getParameterInfo` | `arcpy.Parameter()` | Entire ArcGIS tool parameter system |
| `execute` | `arcpy.AddMessage()` | ArcGIS tool entry-point + logging |
| `check_utm_projection` | `arcpy.Describe()`, `arcpy.AddWarning()`, `arcpy.AddMessage()` | CRS inspection (e.g. `pyproj`, GDAL) + logging |
| `read_strike_dip_data` | `arcpy.Describe()`, `arcpy.da.SearchCursor()` | Feature class I/O (e.g. `geopandas.read_file()`) |
| `extract_z_from_dem` | `arcpy.PointGeometry()`, `arcpy.Point()`, `arcpy.management.GetCellValue()` | Point-at-raster sampling (e.g. `rasterio.sample()`) |
| `extract_dem_cell_centroids` | `arcpy.Describe()`, `arcpy.sa.ExtractByMask()`, `arcpy.RasterToNumPyArray()`, `arcpy.Delete_management()` | Raster clipping + numpy export (e.g. `rasterio` + `shapely` mask) |
| `create_wedge_assignment_raster` | `arcpy.Describe()`, `arcpy.sa.ExtractByMask()`, `arcpy.RasterToNumPyArray()`, `arcpy.NumPyArrayToRaster()`, `arcpy.management.DefineProjection()`, `arcpy.Delete_management()` | Raster I/O (e.g. `rasterio.open()` write) |
| `create_strat_height_raster` | Same as above | Same as above |
| `create_output_strike_dip_fc` | `arcpy.Describe()`, `arcpy.management.CreateFeatureclass()`, `arcpy.ListFields()`, `arcpy.management.AddField()`, `arcpy.da.SearchCursor()`, `arcpy.da.InsertCursor()`, `arcpy.da.UpdateCursor()` | Feature class creation + field I/O (e.g. `geopandas` + `fiona`) |

### Partially arcpy-dependent — core logic is portable

These functions use arcpy **only** for progress reporting and log messages (`arcpy.AddMessage`, `arcpy.AddWarning`, `arcpy.SetProgressor`). All computation is pure numpy. Replacing those calls with `print()` or a standard logging handler makes them fully environment-agnostic.

| Function | arcpy used for | Portable core |
|---|---|---|
| `assign_dem_cells_to_wedges` | Progress bar + messages | Vectorized numpy containment tests |
| `calculate_dem_stratigraphic_heights` | Messages + warnings + progress | Height calculation from wedge geometry |
| `calculate_stratigraphic_heights_in_profile` | Messages + warnings | Height propagation through measurement chain |
| `create_downplunge_view` | Messages + warnings (also calls arcpy-dependent methods) | matplotlib rendering pipeline |

### Fully portable — no arcpy dependency

These functions are pure Python, numpy, or matplotlib and can be used in any environment without modification.

| Function | Dependencies | Notes |
|---|---|---|
| `calculate_fold_axis` | `numpy` | Eigenvalue decomposition on bedding poles |
| `calculate_projected_attitudes` | `numpy` | 3D→2D profile projection and attitude vectors |
| `analyze_wedges` | `numpy` | Parallel vs. intersecting wedge classification |
| `find_line_intersection` | `math` (built-in) | 2D line intersection via Cramer's rule |
| `find_line_intersection_with_params` | `math` (built-in) | Same, with t₁/t₂ parameter output |
| `add_scale_bar` | `matplotlib` | Scale bar renderer for the downplunge view plot |

### Portability summary

The **geological algorithms** (fold axis calculation, profile projection, wedge analysis, height propagation, line intersection) are entirely independent of ArcGIS and form a self-contained computation layer. The **I/O layer** (reading feature classes and rasters, writing outputs, tool registration) is fully arcpy-coupled. A port to another GIS environment would require rewriting only the I/O layer while leaving the geometry and math functions untouched.

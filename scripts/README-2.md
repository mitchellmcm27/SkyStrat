# ThreePointProblemToolbox.pyt — Three-Point Problem Toolbox

An ArcGIS Python Toolbox (`.pyt`) providing two structural geology tools:

1. **Compute Strike/Dip from Point Triplets** — derives a planar orientation (strike, dip) from three points of known 3D position on a bedding surface.
2. **Project Plane Across DEM (Outcrop Trace)** — given a point with a known strike/dip, traces the theoretical line where that plane intersects the topography.

---

## Table of Contents

1. [Overview](#overview)
2. [Dependencies](#dependencies)
3. [Tool 1: Compute Strike/Dip from Point Triplets](#tool-1-compute-strikedip-from-point-triplets)
   - [Input Parameters](#input-parameters)
   - [Output Parameters](#output-parameters)
   - [Workflow](#workflow)
   - [Duplicate Detection](#duplicate-detection)
4. [Tool 2: Project Plane Across DEM (Outcrop Trace)](#tool-2-project-plane-across-dem-outcrop-trace)
   - [Input Parameters](#input-parameters-1)
   - [Output Parameters](#output-parameters-1)
   - [Workflow](#workflow-1)
5. [Key Algorithms](#key-algorithms)
   - [Strike/Dip from Three Points](#1-strikedip-from-three-points-cross-product-method)
   - [Outcrop Trace via Contour](#2-outcrop-trace-via-difference-raster-and-contour)
6. [Helper Functions Reference](#helper-functions-reference)

---

## Overview

The **three-point problem** is a classical structural geology technique: given three points lying on the same geological plane (e.g., the top of a formation) with known X, Y, and Z positions, the orientation of that plane (strike and dip) can be calculated uniquely. This toolbox automates that calculation within ArcGIS Pro, drawing elevations from a DEM when Z values are not already embedded in the features.

The companion tool takes the computed (or user-supplied) orientation and computes the **outcrop trace** — the line where that infinite plane intersects the actual land surface — within a user-defined radius.

---

## Dependencies

| Library | Purpose |
|---|---|
| `arcpy` | ArcGIS spatial operations, raster/feature I/O |
| `numpy` | Plane raster construction for outcrop trace |
| `math`, `hashlib`, `os` | Geometry math, duplicate hashing, file paths |
| `contextlib` | `@contextmanager` for scratch workspace management |

Requires **ArcGIS Pro 3.x**. Spatial Analyst is optional for Tool 1 (speeds up elevation sampling) but **required** for Tool 2 (uses `ExtractByMask` and `Contour`).

---

## Tool 1: Compute Strike/Dip from Point Triplets

**Label**: `Compute Strike/Dip from Point Triplets`
**Class**: `ThreePointProblem`
**Description**: Groups input points into triplets, samples elevation from a DEM, computes the plane normal via cross product, and writes one output point per triplet with the calculated strike and dip.

### Input Parameters

| # | Name | Type | Required | Description |
|---|---|---|---|---|
| 0 | Point Triplets Layer | GPFeatureLayer (Point) | Yes | Input point features forming the triplets |
| 1 | Output Strike/Dip Layer | GPFeatureLayer (Point) | Yes | Existing output point layer to append results into |
| 2 | Elevation Raster | DERasterDataset / GPRasterLayer | Yes | DEM used to sample Z values for each input point |
| 3 | Triplet ID Field | Field | No | Field grouping points into triplets; if omitted, points are grouped sequentially in OID order |

### Output Parameters

The output is appended to the existing **Output Strike/Dip Layer** (parameter 1). The following fields are added if not already present:

| Field | Type | Description |
|---|---|---|
| `STRIKE` | Double | Calculated strike (0–360°, RHR convention) |
| `DIP` | Double | Calculated dip (0–90°) |
| `SRC_OID1`, `SRC_OID2`, `SRC_OID3` | Long | OIDs of the three source points (ascending order) |
| `SRC_LAT1/2/3` | Double | WGS84 latitudes of the three source points |
| `SRC_LON1/2/3` | Double | WGS84 longitudes of the three source points |
| `SRC_HASH` | Text(64) | SHA-256 hash of the triplet's WGS84 coordinates used for deduplication |

Each output point is placed at the **XY position of the second (middle OID) point** in the triplet.

### Workflow

```
1. Resolve input layer paths via arcpy.Describe
2. Validate geometry types (both layers must be Point)
3. Ensure output fields exist; read existing SRC_HASH values to build deduplication set
4. Group input points into triplets:
   a. By Triplet ID Field (all groups must have exactly 3 points), or
   b. By sequential OID order in batches of 3 (total must be divisible by 3)
5. Sample DEM elevations for all input points:
   a. If Spatial Analyst available → ExtractValuesToPoints (interpolated)
   b. Otherwise → numpy cell-index lookup from Raster object
6. For each triplet:
   a. Build (X, Y, Z) tuples from XY geometry + sampled elevation
   b. Project each point to WGS84 and compute stable hash
   c. Skip triplet if hash already exists in output (deduplication)
   d. Compute strike/dip via cross product of plane vectors
   e. Warn and skip if triplet is degenerate (collinear/coincident points)
   f. Insert one output point with strike, dip, source OIDs, lat/lons, and hash
7. Report inserted and skipped counts
```

### Duplicate Detection

Each triplet is identified by a **SHA-256 hash** of its three WGS84 (lat, lon) coordinate pairs. Coordinates are:
- Rounded to 6 decimal places (~0.1 m precision)
- Sorted before hashing (order-independent)

This allows the tool to be run incrementally: new points added to the input layer will produce new results; previously processed triplets are silently skipped.

---

## Tool 2: Project Plane Across DEM (Outcrop Trace)

**Label**: `Project Plane Across DEM (Outcrop Trace)`
**Class**: `ProjectPlaneAcrossDEM`
**Description**: For each input point with a known strike/dip, computes the intersection of the geological plane with the DEM surface within a specified radius, and appends the resulting polylines to an output feature class.

**Requires Spatial Analyst license** (`ExtractByMask`, `Contour`).

### Input Parameters

| # | Name | Type | Required | Description |
|---|---|---|---|---|
| 0 | Input Points (with Strike/Dip) | GPFeatureLayer (Point) | Yes | Points with strike and dip attributes |
| 1 | Use Selected Records | GPBoolean | No | If `True`, process only selected features (default: `False`) |
| 2 | Strike Field | Field | No | Field with strike values; auto-detected if omitted (looks for `Strike`, `Strike_deg`, `Strk`) |
| 3 | Dip Field | Field | No | Field with dip values; auto-detected if omitted (looks for `Dip`, `Dip_deg`) |
| 4 | DEM Raster | DERasterDataset / GPRasterLayer | Yes | Elevation surface for intersection |
| 5 | Target Resolution | GPDouble | No | Resampling cell size in map units; defaults to native DEM resolution |
| 6 | Radius | GPDouble | Yes | Search radius in map units around each point |
| 7 | Output Polyline Feature Class | DEFeatureClass | Yes | Existing or new polyline FC to append outcrop trace segments |

### Output Parameters

Polyline segments are appended to the **Output Polyline Feature Class** (parameter 7). The following fields are populated:

| Field | Type | Description |
|---|---|---|
| `SRC_POINT_OID` / `SRC_OID` | Long | Original OID of the source strike/dip point |
| `STRIKE` | Double | Strike value from the source point |
| `DIP` | Double | Dip value from the source point |

Field name is `SRC_POINT_OID` for geodatabase outputs and `SRC_OID` for shapefiles (10-char limit).

### Workflow

```
1. Validate input points layer (must be Point geometry)
2. Auto-detect strike/dip fields if not explicitly provided
3. Read DEM cell size and spatial reference
4. Copy working points to scratch GDB, preserving original OID, strike, and dip
5. Reproject working points to DEM spatial reference if needed
6. Ensure output polyline FC exists with required fields
7. Check out Spatial Analyst license
8. For each input point:
   a. Build circular buffer polygon of the specified radius
   b. Clip DEM to the circle (ExtractByMask)
   c. Optionally resample clipped DEM to target resolution (Bilinear)
   d. Sample z0 (elevation at the point center) from clipped DEM
   e. Build a plane raster: z_plane = z0 − tan(dip) × (distance along dip direction)
   f. Compute difference raster: diff = DEM_clip − z_plane
   g. Run Contour(diff, interval=1000000, base_contour=0) to extract the 0-isoline
   h. Clip contour to the circle polygon
   i. Append each resulting polyline segment to output FC with source OID, strike, dip
   j. Release scratch rasters and clear workspace cache
9. Check in Spatial Analyst license
10. Report total segments inserted
```

---

## Key Algorithms

### 1. Strike/Dip from Three Points (Cross-Product Method)

Given three 3D points p1, p2, p3:

```
v1 = p2 - p1                              (edge vector)
v2 = p3 - p1                              (edge vector)
n  = v1 × v2                              (normal to plane)
```

The normal is forced to point upward (`nz > 0`) for consistent orientation.

```
dip_deg     = degrees( atan2(sqrt(nx²+ny²), nz) )
dip_dir_deg = degrees( atan2(nx, ny) )    (azimuth of steepest descent)
strike_deg  = (dip_dir_deg − 90°) mod 360°   (RHR convention)
```

Returns `(NaN, NaN)` for degenerate (collinear or coincident) triplets.

### 2. Outcrop Trace via Difference Raster and Contour

The plane equation anchored at the sampled center elevation z0:

```
dip_dir_rad = radians(strike + 90°)
s           = (X − cx) × sin(dip_dir) + (Y − cy) × cos(dip_dir)   (signed distance along dip direction)
z_plane     = z0 − tan(dip_rad) × s
```

The outcrop trace is the set of locations where `DEM = z_plane`, i.e., where `DEM − z_plane = 0`. This is extracted by running the Spatial Analyst `Contour` tool at base contour 0 on the difference raster.

---

## Helper Functions Reference

| Function | Description |
|---|---|
| `_strike_dip_from_three_points(p1, p2, p3)` | Core algorithm: computes (strike, dip) from three 3D points using cross-product plane fitting |
| `_sample_elevations(points_fc, raster, out_pts)` | Samples DEM elevation for all input points; uses `ExtractValuesToPoints` if Spatial Analyst is available, otherwise falls back to numpy cell indexing |
| `_stable_hash_latlons(lats_lons)` | Computes a stable, order-independent SHA-256 hash of a triplet's WGS84 coordinates for deduplication |
| `_autodetect_fields(fc, strike, dip)` | Finds strike/dip fields by scanning field names against common naming patterns |
| `_ensure_fields(fc, field_specs)` | Adds fields to a feature class if they don't already exist |
| `_ensure_polyline_fc(path, sr)` | Creates a polyline feature class with required fields if the path doesn't exist; returns the canonical SRC OID field name |
| `_make_circle(center_geom, radius)` | Returns a circular buffer polygon around a point geometry |
| `_project_to_wgs84(geom)` | Projects a geometry to WGS84 (EPSG:4326) |
| `_layer_or_path_to_fc(value_as_text)` | Resolves a layer name or path to its ArcGIS catalog path |
| `_create_working_points_fc(wksp, sr)` | Creates a scratch point feature class with SRC_OID, STRIKE, DIP fields |
| `_copy_points_with_src_attrs(src, strike_f, dip_f, dst, oid_f)` | Copies points to a scratch FC, preserving original OID and strike/dip values |
| `_has_spatial_analyst()` | Returns `True` if the Spatial Analyst extension is available |
| `_batched(iterable, n)` | Yields successive chunks of size `n` from an iterable |
| `_temp_workspace()` | Context manager that yields the ArcGIS scratch GDB path |
| `_dipdir_from_strike(strike_deg)` | Returns dip direction = (strike + 90°) mod 360° per RHR convention |
| `_deg2rad(d)` | Converts degrees to radians |

---

## Data Constraints and Assumptions

1. **Projected coordinate system recommended**: Elevation sampling and plane distance calculations assume planar (metric) coordinates. Points and DEM should share the same projected CRS; the tool reprojects working points to the DEM CRS if they differ.
2. **Units in map units**: Radius and elevation are assumed to be in the same linear units as the DEM (typically meters).
3. **Three non-collinear points per triplet**: Collinear or coincident points produce a degenerate (zero-magnitude) normal and are skipped with a warning.
4. **Triplet grouping**: Without a Triplet ID field, the total point count must be exactly divisible by 3 and points are processed in ascending OID order.
5. **Spatial Analyst required for Tool 2**: `ExtractByMask` and `Contour` are Spatial Analyst operations; the license is checked out for the duration of each point's processing.
6. **Append semantics**: Both tools **append** to existing output layers rather than overwriting them. For Tool 1, deduplication via `SRC_HASH` prevents double-counting; for Tool 2, all segments are appended unconditionally.

---

## ArcGIS Dependency Analysis

The table below classifies every function and class by how tightly it is coupled to `arcpy`. Functions with no arcpy dependency can be extracted and reused in other GIS environments (QGIS/PyQGIS, GDAL/OGR, GeoPandas, standalone scripts, etc.) with little or no modification.

### Fully arcpy-dependent — ArcGIS required

These functions rely on arcpy for data I/O, geometry operations, or tool-interface registration. They cannot run outside ArcGIS without replacing those arcpy calls with equivalents from another library.

| Function / Class | arcpy APIs used | What would need replacing |
|---|---|---|
| `Toolbox` | `arcpy` toolbox class interface | Entire ArcGIS tool registration system |
| `ThreePointProblem.getParameterInfo` | `arcpy.Parameter()` | Tool parameter system |
| `ThreePointProblem.execute` | `arcpy.Describe()`, `arcpy.da.SearchCursor()`, `arcpy.da.InsertCursor()`, `arcpy.AddMessage()`, `arcpy.AddWarning()`, `arcpy.CreateUniqueName()`, `arcpy.management.GetCount()` | Feature class I/O, logging, scratch name generation |
| `ProjectPlaneAcrossDEM.getParameterInfo` | `arcpy.Parameter()` | Tool parameter system |
| `ProjectPlaneAcrossDEM.execute` | `arcpy.Describe()`, `arcpy.Raster()`, `arcpy.da.SearchCursor/InsertCursor()`, `arcpy.sa.ExtractByMask()`, `arcpy.management.Resample()`, `arcpy.NumPyArrayToRaster()`, `arcpy.management.DefineProjection()`, `arcpy.sa.Contour()`, `arcpy.analysis.Clip()`, `arcpy.GetCellValue_management()`, `arcpy.AddMessage/Warning()`, `arcpy.ClearWorkspaceCache_management()` | Raster clip, resample, contour extraction, I/O — all Spatial Analyst operations |
| `_create_working_points_fc` | `arcpy.CreateUniqueName()`, `arcpy.management.CreateFeatureclass()` | Feature class creation (e.g. `ogr.GetDriverByName` + `CreateDataSource`) |
| `_copy_points_with_src_attrs` | `arcpy.da.SearchCursor()`, `arcpy.da.InsertCursor()` | Feature I/O (e.g. `geopandas.read_file()` + `to_file()`) |
| `_autodetect_fields` | `arcpy.ListFields()`, `arcpy.ExecuteError` | Field introspection (e.g. `layer.fields()` in PyQGIS, `ds.GetLayer().GetLayerDefn()` in OGR) |
| `_make_circle` | arcpy geometry `.buffer()` | Geometry buffering (e.g. `shapely.geometry.Point.buffer()`) |
| `_ensure_polyline_fc` | `arcpy.Exists()`, `arcpy.ListFields()`, `arcpy.management.CreateFeatureclass()`, `arcpy.ExecuteError` | Feature class creation + field inspection |
| `_layer_or_path_to_fc` | `arcpy.Describe()` | Layer-to-path resolution (e.g. `QgsVectorLayer.source()` in QGIS) |
| `_ensure_fields` | `arcpy.ListFields()`, `arcpy.management.AddField()` | Field management (e.g. OGR `FieldDefn`, `CreateField()`) |
| `_project_to_wgs84` | arcpy geometry `.projectAs()`, `arcpy.SpatialReference(4326)` | CRS transformation (e.g. `pyproj.Transformer`, `QgsCoordinateTransform`) |
| `_temp_workspace` | `arcpy.Exists()`, `arcpy.env.scratchGDB`, `arcpy.ExecuteError` | Scratch directory management (e.g. `tempfile.mkdtemp()`) |
| `_has_spatial_analyst` | `arcpy.CheckExtension()` | License check — no equivalent needed outside ArcGIS |
| `_sample_elevations` | `arcpy.sa.ExtractValuesToPoints()`, `arcpy.management.CopyFeatures()`, `arcpy.Raster()`, `arcpy.Describe()`, `arcpy.da.UpdateCursor()`, `arcpy.RasterToNumPyArray()` | Point-at-raster sampling (e.g. `rasterio.sample()`, GDAL `band.ReadAsArray()`) |

### Partially arcpy-dependent — core logic is portable

These contain a self-contained block of pure numpy/math computation embedded inside a larger arcpy-dependent `execute` method. The computation itself is environment-agnostic; only the surrounding I/O and orchestration code is arcpy-specific.

| Code block | arcpy used for | Portable core |
|---|---|---|
| Plane raster construction inside `ProjectPlaneAcrossDEM.execute` (lines ~744–759) | Not used — this block is pure numpy | `np.meshgrid` + trig formula builds the plane elevation grid: `z_plane = z0 − tan(dip) × ((X−cx)·sin(dipdir) + (Y−cy)·cos(dipdir))` |

### Fully portable — no arcpy dependency

These functions are pure Python or pure math and can be used in any environment without modification.

| Function | Dependencies | Notes |
|---|---|---|
| `_strike_dip_from_three_points` | `math` (built-in) | Core algorithm: cross product of edge vectors, derives strike/dip from plane normal |
| `_stable_hash_latlons` | `hashlib` (built-in) | SHA-256 hash of sorted, rounded WGS84 coordinate pairs |
| `_batched` | none | Generator that yields successive chunks of size `n` from an iterable |
| `_deg2rad` | `math` (built-in) | Thin wrapper around `math.radians` |
| `_dipdir_from_strike` | none | RHR arithmetic: `(strike + 90°) mod 360°` |

### Portability summary

The **geological algorithm** — computing strike and dip from three 3D points — is a single, self-contained, arcpy-free function (`_strike_dip_from_three_points`). It depends only on Python's `math` module and can be dropped into any environment as-is. The **plane-intersection math** used for outcrop tracing is also fully portable (pure numpy), but it is currently inlined inside the arcpy-dependent `execute` method rather than extracted as a standalone function.

Everything else — data I/O, geometry operations, raster processing, and tool registration — is tightly coupled to arcpy and the Spatial Analyst extension. A port to another GIS environment (QGIS/PyQGIS, GDAL/rasterio, GeoPandas) would require rewriting the I/O and geometry layers while leaving the two core math blocks untouched.

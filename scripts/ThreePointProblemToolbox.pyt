# -*- coding: utf-8 -*-
# ThreePointProblemToolbox.pyt
# ArcGIS Pro Python Toolbox for computing strike/dip from point triplets on a DEM
# Author: ChatGPT (for Alex T.)
# Tested with ArcGIS Pro 3.x

import arcpy
import os
import math
import hashlib
from contextlib import contextmanager

arcpy.env.overwriteOutput = True

WGS84 = arcpy.SpatialReference(4326)

def _create_working_points_fc(wksp, spatial_ref):
    """Create a scratch points FC with SRC_OID, STRIKE, DIP."""
    name = arcpy.CreateUniqueName("pts_work", wksp)
    arcpy.management.CreateFeatureclass(wksp, os.path.basename(name), "POINT", spatial_reference=spatial_ref)
    _ensure_fields(name, [
        ("SRC_OID", "LONG", None),
        ("STRIKE", "DOUBLE", None),
        ("DIP", "DOUBLE", None),
    ])
    return name

def _copy_points_with_src_attrs(src_layer_or_fc, strike_field, dip_field, dst_fc, oid_field):
    """Copy features from src into dst_fc, writing original OID into SRC_OID and strike/dip into fields."""
    with arcpy.da.SearchCursor(src_layer_or_fc, [oid_field, "SHAPE@", strike_field, dip_field]) as scur, \
         arcpy.da.InsertCursor(dst_fc, ["SHAPE@", "SRC_OID", "STRIKE", "DIP"]) as icur:
        for oid, shp, strike, dip in scur:
            if strike is None or dip is None:
                continue
            icur.insertRow([shp, int(oid), float(strike), float(dip)])

def _autodetect_fields(fc, strike_field_in, dip_field_in):
    flds = {f.name.lower(): f.name for f in arcpy.ListFields(fc)}
    strike = strike_field_in or flds.get("strike") or flds.get("strike_deg") or flds.get("strk") or None
    dip = dip_field_in or flds.get("dip") or flds.get("dip_deg") or None
    if not strike or not dip:
        raise arcpy.ExecuteError("Strike/Dip fields not provided and could not be auto-detected (looked for 'Strike'/'Dip').")
    return strike, dip

def _deg2rad(d):
    return math.radians(float(d))

def _dipdir_from_strike(strike_deg):
    # RHR convention: dip direction = strike + 90
    return (float(strike_deg) + 90.0) % 360.0

def _make_circle(center_geom, radius):
    # Returns a circular polygon (geodesic not required in projected maps; assume planar)
    sr = center_geom.spatialReference
    return center_geom.buffer(radius)

def _ensure_polyline_fc(path, spatial_ref):
    """
    Ensure the polyline featureclass exists with fields for SRC_OID/STRIKE/DIP.
    Returns a tuple (out_fc_path, src_oid_field_name) where src_oid_field_name
    is the actual field to write the source OID into (respects shapefile 10-char limit).
    """
    # If path already exists, find an appropriate SRC OID field (or add one)
    if arcpy.Exists(path):
        names = {f.name.upper(): f.name for f in arcpy.ListFields(path)}
        for cand in ("SRC_POINT_OID", "SRC_OID", "SRCID", "SOURCE_OID"):
            if cand in names:
                return (path, names[cand])
        # add appropriate field name if not present
        # decide short name for shapefile vs long name for gdb
        # detect shapefile by path ending or workspace ending with .gdb
        ws, name = os.path.split(path)
        is_shp = name.lower().endswith(".shp") or (not ws.lower().endswith(".gdb") and os.path.isdir(ws))
        src_field = "SRC_OID" if is_shp else "SRC_POINT_OID"
        _ensure_fields(path, [(src_field, "LONG", None), ("STRIKE", "DOUBLE", None), ("DIP", "DOUBLE", None)])
        return (path, src_field)

    # If missing, create in workspace
    ws, name = os.path.split(path)
    if not ws:
        raise arcpy.ExecuteError(f"Output path '{path}' must include a workspace (e.g., 'C:/folder/out.shp' or 'C:/path/my.gdb/out_fc').")

    # Decide whether to create a shapefile (in a folder) or a geodatabase featureclass
    name_lower = name.lower()
    ws_lower = ws.lower()
    if name_lower.endswith(".shp"):
        # explicit shapefile path provided
        base = os.path.splitext(name)[0]
        arcpy.management.CreateFeatureclass(ws, base + ".shp", "POLYLINE", spatial_reference=spatial_ref)
        out_fc = os.path.join(ws, base + ".shp")
        src_field = "SRC_OID"
    elif ws_lower.endswith(".gdb"):
        # target is a file geodatabase
        arcpy.management.CreateFeatureclass(ws, name, "POLYLINE", spatial_reference=spatial_ref)
        out_fc = os.path.join(ws, name)
        src_field = "SRC_POINT_OID"
    else:
        # workspace isn't a .gdb path; assume folder shapefile
        # create shapefile with provided name (add .shp if missing)
        base = os.path.splitext(name)[0]
        arcpy.management.CreateFeatureclass(ws, base + ".shp", "POLYLINE", spatial_reference=spatial_ref)
        out_fc = os.path.join(ws, base + ".shp")
        src_field = "SRC_OID"

    # Ensure the expected fields exist on the created FC
    _ensure_fields(out_fc, [(src_field, "LONG", None), ("STRIKE", "DOUBLE", None), ("DIP", "DOUBLE", None)])
    return (out_fc, src_field)


def _layer_or_path_to_fc(value_as_text):
    """
    Accepts a GPFeatureLayer parameter's valueAsText (layer or path) and returns a catalog path.
    Works for direct feature class paths too.
    """
    try:
        desc = arcpy.Describe(value_as_text)
        # For layers, Describe.catalogPath gives the source FC; for FC paths, it returns the same path.
        return desc.catalogPath
    except Exception:
        return value_as_text

def _ensure_fields(fc, fields_specs):
    """Add fields if they don't exist. fields_specs = [(name, type, length), ...]"""
    existing = {f.name.upper(): f for f in arcpy.ListFields(fc)}
    to_add = []
    for name, ftype, flen in fields_specs:
        if name.upper() not in existing:
            if ftype.upper() == "TEXT":
                to_add.append((name, ftype, flen))
            else:
                to_add.append((name, ftype, None))
    for name, ftype, flen in to_add:
        if ftype.upper() == "TEXT":
            arcpy.management.AddField(fc, name, ftype, field_length=flen)
        else:
            arcpy.management.AddField(fc, name, ftype)

def _project_to_wgs84(geom):
    # geom is a Geometry (PointGeometry)
    sr = getattr(geom, "spatialReference", None)
    if sr and sr.factoryCode == 4326:
        return geom
    return geom.projectAs(WGS84)

def _stable_hash_latlons(lats_lons):
    """lats_lons: list of (lat, lon) for the 3 points. Order-independent hash (sort first)."""
    # round to ~1e-6 deg (~0.1 m) to make hashing robust to tiny reprojection differences
    rounded = sorted([(round(lat, 6), round(lon, 6)) for (lat, lon) in lats_lons])
    s = ";".join([f"{lat:.6f},{lon:.6f}" for lat, lon in rounded])
    return hashlib.sha256(s.encode("utf-8")).hexdigest()[:64]  # 64 hex chars

def _strike_dip_from_three_points(p1, p2, p3):
    """
    p1,p2,p3: (x,y,z). Returns (strike_deg_0_360, dip_deg_0_90).
    Order-invariant: forces nz > 0, then derives down-dip and strike by RHR.
    """
    # 1) vectors in 3D
    v1 = (p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2])
    v2 = (p3[0]-p1[0], p3[1]-p1[1], p3[2]-p1[2])

    # 2) unnormalized normal n = v1 x v2
    nx = v1[1]*v2[2] - v1[2]*v2[1]
    ny = v1[2]*v2[0] - v1[0]*v2[2]
    nz = v1[0]*v2[1] - v1[1]*v2[0]

    # Degenerate?
    norm = math.sqrt(nx*nx + ny*ny + nz*nz)
    if norm == 0:
        return (float("nan"), float("nan"))

    # 3) Force consistent orientation: nz > 0 (upward)
    if nz < 0:
        nx, ny, nz = -nx, -ny, -nz

    # 4) Dip: angle from horizontal = atan2(|horizontal n|, vertical n)
    h = math.sqrt(nx*nx + ny*ny)           # horizontal component magnitude
    dip_deg = math.degrees(math.atan2(h, nz))  # nz >= 0 here

    # 5) Down-dip azimuth from horizontal projection of the normal
    # For nz>0, the direction of steepest *descent* is (nx, ny)
    dipdir = math.degrees(math.atan2(nx, ny))  # atan2(east, north)
    if dipdir < 0:
        dipdir += 360.0

    # 6) Strike (RHR): 90° counterclockwise from dip direction
    strike = (dipdir - 90.0) % 360.0

    return (strike, dip_deg)

def _batched(iterable, n):
    batch = []
    for item in iterable:
        batch.append(item)
        if len(batch) == n:
            yield batch
            batch = []
    if batch:
        yield batch

@contextmanager
def _temp_workspace():
    if not arcpy.Exists(arcpy.env.scratchGDB):
        raise arcpy.ExecuteError("Scratch GDB not available. Check your ArcGIS Pro environment.")
    yield arcpy.env.scratchGDB

def _has_spatial_analyst():
    try:
        return arcpy.CheckExtension("Spatial") in ("Available", "Enabled")
    except:
        return False

def _sample_elevations(points_fc, raster, out_pts):
    """
    Adds elevation value to each point feature in points_fc and writes to out_pts.
    Prefers Spatial Analyst's ExtractValuesToPoints; otherwise does numpy sampler.
    Returns the name of the numeric field containing elevation (e.g., 'RASTERVALU' or 'ELEV').
    """
    if _has_spatial_analyst():
        # Use Spatial Analyst if available
        try:
            arcpy.CheckOutExtension("Spatial")
            arcpy.sa.ExtractValuesToPoints(
                points_fc,
                raster,
                out_pts,
                interpolate_values="INTERPOLATE",
                add_attributes="VALUE_ONLY"
            )
        finally:
            # Always give the license back
            arcpy.CheckInExtension("Spatial")
        return "RASTERVALU"
    else:
        # --- existing license-free fallback code stays the same below ---
        arcpy.management.CopyFeatures(points_fc, out_pts)
        r = arcpy.Raster(raster)
        desc = arcpy.Describe(r)
        cellw = r.meanCellWidth
        cellh = r.meanCellHeight
        extent = desc.Extent
        x0, y0 = extent.XMin, extent.YMax
        band = r
        elev_field = "ELEV"
        _ensure_fields(out_pts, [(elev_field, "DOUBLE", None)])
        with arcpy.da.UpdateCursor(out_pts, ["SHAPE@XY", elev_field]) as cur:
            for (x, y), _ in cur:
                col = int((x - x0) / cellw)
                row = int((y0 - y) / cellh)
                try:
                    val = band.read(1, row, col)
                except:
                    arr = arcpy.RasterToNumPyArray(band, arcpy.Point(x, y), 1, 1)
                    val = float(arr[0, 0])
                cur.updateRow([(x, y), float(val)])
        return elev_field

class Toolbox(object):
    def __init__(self):
        self.label = "Three Point Problem Toolbox"
        self.alias = "threepoint"
        self.tools = [ThreePointProblem, ProjectPlaneAcrossDEM]  # <-- add here

class ThreePointProblem(object):
    def __init__(self):
        self.label = "Compute Strike/Dip from Point Triplets"
        self.description = ("Given a feature class of point triplets and an elevation raster, "
                            "compute strike (0–360) and dip (0–90) for each new triplet and write a point at the XY centroid to the output feature class. "
                            "Triplets already processed are skipped using a WGS84 lat/lon hash stored in SRC_HASH.")
        self.canRunInBackground = True

    def getParameterInfo(self):
        in_pts = arcpy.Parameter(
            displayName="Point Triplets Layer",
            name="in_points",
            datatype="GPFeatureLayer",   # <-- was DEFeatureClass
            parameterType="Required",
            direction="Input"
        )

        out_fc = arcpy.Parameter(
            displayName="Output Strike/Dip Layer",
            name="out_fc",
            datatype="GPFeatureLayer",   # <-- was DEFeatureClass
            parameterType="Required",
            direction="Input"
        )

        in_raster = arcpy.Parameter(
            displayName="Elevation Raster",
            name="in_raster",
            datatype=["DERasterDataset", "GPRasterLayer"],
            parameterType="Required",
            direction="Input"
        )

        triplet_field = arcpy.Parameter(
            displayName="Triplet ID Field (optional)",
            name="triplet_id_field",
            datatype="Field",
            parameterType="Optional",
            direction="Input"
        )
        triplet_field.parameterDependencies = [in_pts.name]

        return [in_pts, out_fc, in_raster, triplet_field]

    def isLicensed(self):
        # Tool runs with or without Spatial Analyst; SA just speeds up sampling.
        return True

    def updateParameters(self, params):
        return

    def updateMessages(self, params):
        return

    def execute(self, params, messages):
        in_points_layer = params[0].valueAsText
        out_fc_layer = params[1].valueAsText
        in_raster = params[2].valueAsText
        triplet_field = params[3].valueAsText if params[3].value else None

        # Resolve GPFeatureLayer selections to catalog paths
        in_points = _layer_or_path_to_fc(in_points_layer)
        out_fc = _layer_or_path_to_fc(out_fc_layer)

        # 0) Ensure output fields exist
        out_fields = [
            ("STRIKE", "DOUBLE", None),
            ("DIP", "DOUBLE", None),
            ("SRC_OID1", "LONG", None),
            ("SRC_OID2", "LONG", None),
            ("SRC_OID3", "LONG", None),
            ("SRC_LAT1", "DOUBLE", None),
            ("SRC_LAT2", "DOUBLE", None),
            ("SRC_LAT3", "DOUBLE", None),
            ("SRC_LON1", "DOUBLE", None),
            ("SRC_LON2", "DOUBLE", None),
            ("SRC_LON3", "DOUBLE", None),
            ("SRC_HASH", "TEXT", 64),
        ]
        _ensure_fields(out_fc, out_fields)

        # 1) Build a set of existing hashes to skip duplicates
        existing_hashes = set()
        with arcpy.da.SearchCursor(out_fc, ["SRC_HASH"]) as cur:
            for (h,) in cur:
                if h:
                    existing_hashes.add(h)

        # 2) Collect input points, grouped into triplets
        oid_field = arcpy.Describe(in_points).OIDFieldName
        sr = arcpy.Describe(in_points).spatialReference

        # Validate geometry type = Point
        in_desc = arcpy.Describe(in_points)
        out_desc = arcpy.Describe(out_fc)
        if getattr(in_desc, "shapeType", "").lower() != "point":
            raise arcpy.ExecuteError("The Point Triplets Layer must be a point feature layer.")
        if getattr(out_desc, "shapeType", "").lower() != "point":
            raise arcpy.ExecuteError("The Output Strike/Dip Layer must be a point feature layer.")

        if triplet_field:
            fields = [oid_field, "SHAPE@", triplet_field]
            triplet_map = {}
            with arcpy.da.SearchCursor(in_points, fields) as cur:
                for oid, geom, tid in cur:
                    if tid is None:
                        arcpy.AddWarning(f"Point OID {oid} has null TripletID; it will be ignored.")
                        continue
                    triplet_map.setdefault(tid, []).append((oid, geom))
            # Validate counts
            bad = [tid for tid, pts in triplet_map.items() if len(pts) != 3]
            if bad:
                raise arcpy.ExecuteError(
                    f"The following TripletID groups do not contain exactly 3 points: {bad[:10]}{'...' if len(bad)>10 else ''}"
                )
            groups = [sorted(pts, key=lambda x: x[0]) for _, pts in triplet_map.items()]
        else:
            # Group by OID order in batches of 3
            fields = [oid_field, "SHAPE@"]
            pts = []
            with arcpy.da.SearchCursor(in_points, fields, sql_clause=(None, f"ORDER BY {oid_field}")) as cur:
                for oid, geom in cur:
                    pts.append((oid, geom))
            if len(pts) % 3 != 0:
                raise arcpy.ExecuteError(
                    f"Input has {len(pts)} points; not divisible by 3. Add/remove points or specify a Triplet ID field."
                )
            groups = list(_batched(pts, 3))

        arcpy.AddMessage(f"Found {len(groups)} triplets in input.")

        # 3) Elevation sampling: write a temp copy with elevation field added
        with _temp_workspace() as scratch_gdb:
            temp_pts = arcpy.CreateUniqueName("pts_with_z", scratch_gdb)
            elev_field_name = _sample_elevations(in_points, in_raster, temp_pts)

            # For fast lookup: OID -> (geom, elev)
            elev_lookup = {}
            fields = [oid_field, "SHAPE@", elev_field_name]
            with arcpy.da.SearchCursor(temp_pts, fields) as cur:
                for oid, geom, z in cur:
                    elev_lookup[oid] = (geom, float(z))

        # 4) Prepare insert cursor into output
        out_fields_write = [
            "SHAPE@XY", "STRIKE", "DIP",
            "SRC_OID1", "SRC_OID2", "SRC_OID3",
            "SRC_LAT1", "SRC_LAT2", "SRC_LAT3",
            "SRC_LON1", "SRC_LON2", "SRC_LON3",
            "SRC_HASH"
        ]

        inserted = 0
        skipped = 0

        with arcpy.da.InsertCursor(out_fc, out_fields_write) as icur:
            for trip in groups:
                # trip: [(oid, geom), ...] length 3
                try:
                    # Pull XY and elevation
                    items = []
                    for oid, _geom in trip:
                        if oid not in elev_lookup:
                            raise RuntimeError(f"Missing elevation for OID {oid}.")
                        g, z = elev_lookup[oid]
                        items.append((oid, g, z))

                    # WGS84 lat/lon for hashing
                    latlons = []
                    for oid, g, _z in items:
                        gw_geom = _project_to_wgs84(g)    # <- pass geometry
                        pt = gw_geom.centroid             # Point (has X/Y)
                        latlons.append((pt.Y, pt.X))
                    src_hash = _stable_hash_latlons(latlons)

                    if src_hash in existing_hashes:
                        skipped += 1
                        continue

                    # XY centroid (in source SR)
                    xs = [g.centroid.X for _, g, _ in items]
                    ys = [g.centroid.Y for _, g, _ in items]
                    #cx, cy = sum(xs) / 3.0, sum(ys) / 3.0
                    cx = xs[1]
                    cy = ys[1]

                    # 3D points (use original XY + sampled Z)
                    p1 = (items[0][1].centroid.X, items[0][1].centroid.Y, items[0][2])
                    p2 = (items[1][1].centroid.X, items[1][1].centroid.Y, items[1][2])
                    p3 = (items[2][1].centroid.X, items[2][1].centroid.Y, items[2][2])

                    strike, dip = _strike_dip_from_three_points(p1, p2, p3)
                    if math.isnan(strike) or math.isnan(dip):
                        arcpy.AddWarning(
                            f"Degenerate triplet (collinear or identical points); skipping OIDs {[it[0] for it in items]}."
                        )
                        continue

                    # Order lat/lon arrays to match SRC_OID1..3 ascending by OID
                    items_sorted = sorted(items, key=lambda x: x[0])
                    oids_sorted = [it[0] for it in items_sorted]

                    lat_map = {}
                    lon_map = {}
                    for idx, (oid, g, _z) in enumerate(items_sorted, start=1):
                        gw_geom = _project_to_wgs84(g)    # <- pass geometry
                        pt = gw_geom.centroid
                        lat_map[idx] = pt.Y
                        lon_map[idx] = pt.X

                    row = [
                        (cx, cy), strike, dip,
                        oids_sorted[0], oids_sorted[1], oids_sorted[2],
                        lat_map[1], lat_map[2], lat_map[3],
                        lon_map[1], lon_map[2], lon_map[3],
                        src_hash
                    ]
                    icur.insertRow(row)
                    existing_hashes.add(src_hash)
                    inserted += 1
                except Exception as ex:
                    arcpy.AddWarning(f"Failed triplet {[oid for oid, _ in trip]}: {ex}")

        arcpy.AddMessage(f"Inserted {inserted} new strike/dip points. Skipped {skipped} already-processed triplets.")

class ProjectPlaneAcrossDEM(object):
    def __init__(self):
        self.label = "Project Plane Across DEM (Outcrop Trace)"
        self.description = ("Given a point layer with strike/dip and a DEM, computes the polyline where the plane "
                            "intersects the topography within a radius around each point. Appends to an output line FC.")
        self.canRunInBackground = True

    def getParameterInfo(self):
        pts = arcpy.Parameter(
            displayName="Input Points (with Strike/Dip)",
            name="in_points",
            datatype="GPFeatureLayer",   # dropdown of map layers
            parameterType="Required",
            direction="Input"
        )

        use_sel = arcpy.Parameter(
            displayName="Use selected records",
            name="use_selected",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input"
        )
        use_sel.value = False

        strike_f = arcpy.Parameter(
            displayName="Strike Field (optional)",
            name="strike_field",
            datatype="Field",
            parameterType="Optional",
            direction="Input"
        )
        strike_f.parameterDependencies = [pts.name]

        dip_f = arcpy.Parameter(
            displayName="Dip Field (optional)",
            name="dip_field",
            datatype="Field",
            parameterType="Optional",
            direction="Input"
        )
        dip_f.parameterDependencies = [pts.name]

        dem = arcpy.Parameter(
            displayName="DEM Raster",
            name="in_dem",
            datatype=["DERasterDataset", "GPRasterLayer"],
            parameterType="Required",
            direction="Input"
        )

        res = arcpy.Parameter(
            displayName="Target Resolution (map units; optional)",
            name="resolution",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input"
        )

        rad = arcpy.Parameter(
            displayName="Radius (map units)",
            name="radius",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input"
        )

        out_fc = arcpy.Parameter(
            displayName="Output Polyline Feature Class (append if exists)",
            name="out_fc",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Input"   # <-- was "Output"; with "Input" we won't auto-delete
        )

        return [pts, use_sel, strike_f, dip_f, dem, res, rad, out_fc]

    def isLicensed(self):
        return True

    def execute(self, params, messages):
        """
        Execute method for ProjectPlaneAcrossDEM using Spatial Analyst Contour.
        Param order expected:
        0: in_points (GPFeatureLayer / path)
        1: use_selected (GPBoolean)
        2: strike_field (Field) optional
        3: dip_field (Field) optional
        4: in_dem (DERasterDataset / GPRasterLayer)
        5: resolution (GPDouble) optional
        6: radius (GPDouble)
        7: out_fc (DEFeatureClass path) - append target
        """
        import time
        import numpy as np
        import gc

        # --- Read params ---
        in_points_layer = params[0].valueAsText
        use_selected = bool(params[1].value) if params[1].value is not None else False
        strike_field = params[2].valueAsText if params[2].value else None
        dip_field = params[3].valueAsText if params[3].value else None
        in_dem = params[4].valueAsText
        resolution = float(params[5].valueAsText) if params[5].value else None
        radius = float(params[6].valueAsText)
        out_fc_param = params[7].valueAsText

        arcpy.AddMessage("[Start] ProjectPlaneAcrossDEM (Contour) starting...")
        start_all = time.time()

        # Resolve layer path if needed (keep selection behavior)
        try:
            desc_pts_layer = arcpy.Describe(in_points_layer)
        except Exception as e:
            raise arcpy.ExecuteError(f"Cannot access input points layer: {e}")
        if getattr(desc_pts_layer, "shapeType", "").lower() != "point":
            raise arcpy.ExecuteError("Input points must be a point layer/feature class.")

        # Autodetect strike/dip fields if not supplied
        strike_field, dip_field = _autodetect_fields(in_points_layer, strike_field, dip_field)

        # DEM info
        dem_r = arcpy.Raster(in_dem)
        dem_sr = dem_r.spatialReference
        dem_cx = dem_r.meanCellWidth
        dem_cy = dem_r.meanCellHeight
        base_cs = max(dem_cx, dem_cy)
        target_cs = resolution if (resolution and resolution > 0) else base_cs

        arcpy.AddMessage(f"[Setup] radius={radius}, requested resolution={resolution}, using cellsize={target_cs}")

        # Materialize working points (preserve original OID in SRC_OID)
        oid_field = desc_pts_layer.OIDFieldName
        in_sr = desc_pts_layer.spatialReference

        with _temp_workspace() as wksp:
            if use_selected:
                arcpy.AddMessage("[Setup] Copying selected features from layer to scratch...")
                sel_fc = arcpy.CreateUniqueName("pts_sel", wksp)
                arcpy.management.CopyFeatures(in_points_layer, sel_fc)
                if int(arcpy.management.GetCount(sel_fc).getOutput(0)) == 0:
                    arcpy.AddWarning("No selected records found; exiting.")
                    return
                src_for_copy = sel_fc
            else:
                arcpy.AddMessage("[Setup] Materializing all input points to scratch...")
                src_for_copy = in_points_layer

            # Create working points FC that stores SRC_OID, STRIKE, DIP
            pts_work = arcpy.CreateUniqueName("pts_work", wksp)
            arcpy.management.CreateFeatureclass(wksp, os.path.basename(pts_work), "POINT", spatial_reference=in_sr)
            _ensure_fields(pts_work, [("SRC_OID", "LONG", None), ("STRIKE", "DOUBLE", None), ("DIP", "DOUBLE", None)])

            copied = 0
            with arcpy.da.SearchCursor(src_for_copy, [oid_field, "SHAPE@", strike_field, dip_field]) as scur, \
                arcpy.da.InsertCursor(pts_work, ["SHAPE@", "SRC_OID", "STRIKE", "DIP"]) as icur_tmp:
                for oid, shp, s, d in scur:
                    if s is None or d is None:
                        arcpy.AddWarning(f"[Setup] Source OID {oid} missing strike/dip; skipping.")
                        continue
                    icur_tmp.insertRow([shp, int(oid), float(s), float(d)])
                    copied += 1
            if copied == 0:
                arcpy.AddWarning("No usable points found (all missing strike/dip). Exiting.")
                return

            # Project working pts into DEM SR if needed
            if in_sr.factoryCode != dem_sr.factoryCode:
                arcpy.AddMessage("[Setup] Projecting working points to DEM spatial reference...")
                pts_fc = arcpy.management.Project(pts_work, arcpy.CreateUniqueName("pts_dem_sr", wksp), dem_sr)
            else:
                pts_fc = arcpy.management.CopyFeatures(pts_work, arcpy.CreateUniqueName("pts_dem_sr", wksp))

        # Ensure output FC exists and get canonical SRC OID field name
        out_fc, src_oid_outfield = _ensure_polyline_fc(out_fc_param, dem_sr)
        arcpy.AddMessage(f"[Setup] Output: {out_fc} (writing original OID into field '{src_oid_outfield}')")

        # Count points
        npts = int(arcpy.management.GetCount(pts_fc).getOutput(0))
        arcpy.AddMessage(f"[Setup] {npts} points to process.")

        # Check Spatial Analyst
        try:
            arcpy.CheckOutExtension("Spatial")
        except Exception as e:
            raise arcpy.ExecuteError(f"Spatial Analyst license required: {e}")

        # Prepare insert cursor for output
        insert_fields = ["SHAPE@", src_oid_outfield, "STRIKE", "DIP"]
        inserted_total = 0

        # Open insert cursor once
        with arcpy.da.InsertCursor(out_fc, insert_fields) as icur:
            pcounter = 0
            # iterate through working points
            with arcpy.da.SearchCursor(pts_fc, ["SRC_OID", "SHAPE@", "STRIKE", "DIP"]) as cur:
                for src_oid, geom, strike, dip in cur:
                    pcounter += 1
                    arcpy.AddMessage(f"\n[Point {pcounter}/{npts}] SRC_OID={src_oid} strike={strike} dip={dip}")
                    t0 = time.time()

                    try:
                        strike = float(strike); dip = float(dip)
                    except Exception:
                        arcpy.AddWarning(f"[Point {src_oid}] Strike/Dip conversion failed; skipping.")
                        continue

                    # Build circular mask around the point and clip DEM
                    try:
                        circle = _make_circle(geom, radius)
                    except Exception as e:
                        arcpy.AddWarning(f"[Point {src_oid}] Could not build circle buffer: {e}")
                        continue

                    with _temp_workspace() as wksp2:
                        # materialize circle
                        circle_fc = arcpy.CreateUniqueName("circle", wksp2)
                        arcpy.management.CopyFeatures(circle, circle_fc)

                        arcpy.AddMessage(f"[Point {src_oid}] Extracting DEM by mask...")
                        dem_clip_res = arcpy.sa.ExtractByMask(dem_r, circle_fc)
                        dem_clip_path = arcpy.CreateUniqueName("dem_clip", wksp2)
                        dem_clip_res.save(dem_clip_path)
                        dem_clip = arcpy.Raster(dem_clip_path)

                        # Optionally resample
                        if abs(target_cs - base_cs) / base_cs > 0.05:
                            arcpy.AddMessage(f"[Point {src_oid}] Resampling DEM clip to {target_cs} ...")
                            dem_rs_path = arcpy.CreateUniqueName("dem_clip_rs", wksp2)
                            arcpy.management.Resample(dem_clip_path, dem_rs_path, f"{target_cs} {target_cs}", "BILINEAR")
                            dem_clip_path = dem_rs_path
                            dem_clip = arcpy.Raster(dem_clip_path)

                        # Describe clipped DEM and grid geometry
                        dem_desc = arcpy.Describe(dem_clip)
                        ext = dem_desc.extent
                        csx = dem_clip.meanCellWidth
                        csy = dem_clip.meanCellHeight
                        # compute ncols/nrows from extent & cellsize (round to nearest int)
                        ncols = int(round((ext.XMax - ext.XMin) / csx))
                        nrows = int(round((ext.YMax - ext.YMin) / csy))

                        arcpy.AddMessage(f"[Point {src_oid}] dem_clip extent = ({ext.XMin:.3f}, {ext.YMin:.3f}) - ({ext.XMax:.3f}, {ext.YMax:.3f}), cols={ncols}, rows={nrows}, cs=({csx:.3f},{csy:.3f})")

                        # Sample z0 at the point location
                        center = geom.centroid if hasattr(geom, "centroid") else geom
                        try:
                            z0 = float(arcpy.GetCellValue_management(dem_clip, f"{center.X} {center.Y}").getOutput(0))
                        except Exception as e:
                            arcpy.AddWarning(f"[Point {src_oid}] Could not sample DEM at center: {e}; skipping.")
                            continue
                        arcpy.AddMessage(f"[Point {src_oid}] z0 = {z0:.3f}")

                        # Build plane raster values on same grid as dem_clip
                        x_left = ext.XMin
                        y_top = ext.YMax
                        cols = np.arange(ncols, dtype=float)
                        rows = np.arange(nrows, dtype=float)
                        Xc = x_left + (cols + 0.5) * csx
                        Yc = y_top  - (rows + 0.5) * csy
                        Xg, Yg = np.meshgrid(Xc, Yc)  # shape (nrows, ncols)

                        dipdir = (strike + 90.0) % 360.0
                        diprad = math.radians(dip)
                        dd_rad = math.radians(dipdir)
                        ux = math.sin(dd_rad)
                        uy = math.cos(dd_rad)
                        slope = math.tan(diprad)

                        # signed distance along dip direction from the center
                        s = (Xg - center.X) * ux + (Yg - center.Y) * uy
                        z_plane = z0 - slope * s

                        # Convert plane numpy array to raster with the same georeferencing
                        plane_arr = np.array(z_plane, dtype=float)
                        origin_pt = arcpy.Point(ext.XMin, ext.YMin)
                        plane_raster = arcpy.NumPyArrayToRaster(plane_arr, origin_pt, csx, csy, value_to_nodata=np.nan)
                        plane_path = arcpy.CreateUniqueName("plane", wksp2)
                        plane_raster.save(plane_path)
                        # define projection to DEM SR
                        try:
                            arcpy.management.DefineProjection(plane_path, dem_sr)
                        except Exception:
                            pass

                        # Compute diff raster using Spatial Analyst arithmetic
                        arcpy.AddMessage(f"[Point {src_oid}] Computing diff raster (DEM_clip - plane)...")
                        dem_r_obj = arcpy.Raster(dem_clip_path)
                        plane_r_obj = arcpy.Raster(plane_path)
                        diff_r_obj = dem_r_obj - plane_r_obj
                        diff_path = arcpy.CreateUniqueName("diff_r", wksp2)
                        diff_r_obj.save(diff_path)
                        try:
                            arcpy.management.DefineProjection(diff_path, dem_sr)
                        except Exception:
                            pass
                        arcpy.AddMessage(f"[Point {src_oid}] diff raster saved -> {diff_path}")

                        # Run Spatial Analyst Contour to produce 0-contour polyline
                        contour_fc = None
                        try:
                            arcpy.AddMessage(f"[Point {src_oid}] Running Spatial Analyst Contour (0.0) ...")
                            temp_contour = arcpy.CreateUniqueName("contour_raw", wksp2)
                            # contour_interval set to 1 (we extract base_contour=0)
                            arcpy.sa.Contour(arcpy.Raster(diff_path), temp_contour, 1000000, 0)
                            contour_fc = temp_contour
                            arcpy.AddMessage(f"[Point {src_oid}] Contour succeeded -> {contour_fc}")
                        except Exception as ce:
                            arcpy.AddWarning(f"[Point {src_oid}] Contour failed: {ce}. No contour will be produced for this point.")
                            contour_fc = None

                        # If a contour exists, clip to circle and append
                        if contour_fc and arcpy.Exists(contour_fc):
                            clipped = arcpy.CreateUniqueName("trace_clip", wksp2)
                            try:
                                arcpy.analysis.Clip(contour_fc, circle_fc, clipped)
                            except Exception as ce:
                                arcpy.AddWarning(f"[Point {src_oid}] Clip failed: {ce}. Attempting to append contour directly.")
                                clipped = contour_fc

                            cnt_lines = int(arcpy.management.GetCount(clipped).getOutput(0))
                            arcpy.AddMessage(f"[Point {src_oid}] clipped trace features: {cnt_lines}")
                            if cnt_lines == 0:
                                arcpy.AddWarning(f"[Point {src_oid}] No contour segments after clipping. Try increasing radius or check DEM coverage.")
                            else:
                                added = 0
                                with arcpy.da.SearchCursor(clipped, ["SHAPE@"]) as lcur:
                                    for (ln,) in lcur:
                                        icur.insertRow([ln, int(src_oid), float(strike), float(dip)])
                                        added += 1
                                inserted_total += added
                                arcpy.AddMessage(f"[Point {src_oid}] appended {added} segments.")
                        else:
                            arcpy.AddWarning(f"[Point {src_oid}] No contour produced for this point.")

                        # Cleanup objects & avoid locks
                        try:
                            del dem_r_obj, plane_r_obj, diff_r_obj
                        except:
                            pass
                        try:
                            gc.collect()
                            arcpy.ClearWorkspaceCache_management()
                        except:
                            pass

                        arcpy.AddMessage(f"[Point {src_oid}] elapsed {time.time() - t0:.2f}s")

                # end for points

        # Check in Spatial Analyst
        try:
            arcpy.CheckInExtension("Spatial")
        except:
            pass

        arcpy.AddMessage(f"[Finish] Inserted {inserted_total} segments across {npts} points. Total elapsed: {time.time() - start_all:.2f}s")
        return

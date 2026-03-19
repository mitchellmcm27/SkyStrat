# -*- coding: utf-8 -*-
"""
Three Point Problem — QGIS Processing Script
=============================================

Port of ThreePointProblemToolbox.pyt (ArcGIS Pro) to QGIS Processing.
Contains two algorithms:

  1. Compute Strike/Dip from Point Triplets  (ComputeStrikeDipAlgorithm)
  2. Project Plane Across DEM (Outcrop Trace) (ProjectPlaneAcrossDEMAlgorithm)

QGIS discovers all QgsProcessingAlgorithm subclasses in a scripts-folder file,
so both tools appear automatically in the Processing Toolbox.

Installation
------------
Copy this file to your QGIS Processing Scripts folder:
  Windows : %APPDATA%\\QGIS\\QGIS3\\profiles\\default\\processing\\scripts\\
  Linux   : ~/.local/share/QGIS/QGIS3/profiles/default/processing/scripts/
  macOS   : ~/Library/Application Support/QGIS/QGIS3/profiles/default/processing/scripts/

After copying, open the Processing Toolbox and click the Refresh Scripts button.
Both tools appear under: Scripts > Structural Geology

Dependencies: QGIS 3.x, GDAL/OGR (bundled with QGIS), numpy
"""

import gc
import math
import os
import tempfile
import time

import numpy as np
from osgeo import gdal, ogr, osr

from qgis.core import (
    QgsCoordinateReferenceSystem,
    QgsCoordinateTransform,
    QgsFeature,
    QgsFeatureSink,
    QgsField,
    QgsFields,
    QgsGeometry,
    QgsPointXY,
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingException,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterFeatureSink,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterField,
    QgsProcessingParameterNumber,
    QgsProcessingParameterRasterLayer,
    QgsWkbTypes,
)
from qgis.PyQt.QtCore import QVariant

gdal.UseExceptions()


# =============================================================================
# Pure-math helpers — identical to original, no GIS dependency
# =============================================================================

def _deg2rad(d):
    return math.radians(float(d))


def _dipdir_from_strike(strike_deg):
    """RHR convention: dip direction = strike + 90°."""
    return (float(strike_deg) + 90.0) % 360.0

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
    """Yield successive chunks of size n from iterable."""
    batch = []
    for item in iterable:
        batch.append(item)
        if len(batch) == n:
            yield batch
            batch = []
    if batch:
        yield batch


# =============================================================================
# GDAL/QGIS utility helpers — replacements for arcpy I/O functions
# =============================================================================

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

def _autodetect_field(source, provided, candidates):
    """
    Return provided field name if given; otherwise search source.fields() for
    the first match against the candidate names (case-insensitive).
    Replaces _autodetect_fields() from the original.
    """
    if provided:
        return provided
    names_lower = {f.name().lower(): f.name() for f in source.fields()}
    for cand in candidates:
        if cand.lower() in names_lower:
            return names_lower[cand.lower()]
    return None


def _write_circle_shapefile(cx, cy, radius, proj_wkt, path):
    """
    Write a circular buffer polygon around (cx, cy) to a temp shapefile.
    Used as a GDAL cutline for DEM clipping. Replaces arcpy buffer + shapefile export.
    """
    driver = ogr.GetDriverByName('ESRI Shapefile')
    if os.path.exists(path):
        driver.DeleteDataSource(path)
    ds = driver.CreateDataSource(path)
    srs = osr.SpatialReference()
    srs.ImportFromWkt(proj_wkt)
    layer = ds.CreateLayer('circle', srs=srs, geom_type=ogr.wkbPolygon)
    ogr_pt = ogr.Geometry(ogr.wkbPoint)
    ogr_pt.AddPoint(cx, cy)
    circle_geom = ogr_pt.Buffer(radius, 64)    # 64-segment approximation
    feat = ogr.Feature(layer.GetLayerDefn())
    feat.SetGeometry(circle_geom)
    layer.CreateFeature(feat)
    ds.FlushCache()
    ds = None
    return circle_geom   # return OGR geometry for later intersection use


# =============================================================================
# Algorithm 1: Compute Strike/Dip from Point Triplets
# =============================================================================

class ComputeStrikeDipAlgorithm(QgsProcessingAlgorithm):
    """
    Given a point feature class of triplets and an elevation DEM, computes
    strike (0-360°, RHR) and dip (0-90°) for each triplet and writes one
    output point per triplet.

    Direct port of the ThreePointProblem tool from ThreePointProblemToolbox.pyt.
    """

    IN_POINTS     = 'in_points'
    OUT_FC        = 'out_fc'
    IN_RASTER     = 'in_raster'
    TRIPLET_FIELD = 'triplet_id_field'

    # -------------------------------------------------------------------------
    # Algorithm identity
    # -------------------------------------------------------------------------

    def createInstance(self):
        return ComputeStrikeDipAlgorithm()

    def name(self):
        return 'computestrikedip'

    def displayName(self):
        return 'Compute Strike/Dip from Point Triplets'

    def group(self):
        return 'Structural Geology'

    def groupId(self):
        return 'structuralgeology'

    def shortHelpString(self):
        return (
            'Given a point layer of triplets and an elevation DEM, computes '
            'strike (0–360°, RHR convention) and dip (0–90°) for each triplet '
            'by fitting a plane to three 3D points via cross product.\n\n'
            'One output point is placed at the location of the middle-OID point '
            'in each triplet.\n\n'
            'Grouping options:\n'
            '  • Triplet ID field — every group must have exactly 3 points.\n'
            '  • Sequential order — total feature count must be divisible by 3.'
        )

    # -------------------------------------------------------------------------
    # Parameter definitions — mirrors getParameterInfo() in the arcpy version
    # -------------------------------------------------------------------------

    def initAlgorithm(self, config=None):
        # 0 – Input point triplets
        self.addParameter(QgsProcessingParameterFeatureSource(
            self.IN_POINTS,
            'Point Triplets Layer',
            [QgsProcessing.TypeVectorPoint],
        ))

        # 1 – Output strike/dip points (new layer each run)
        self.addParameter(QgsProcessingParameterFeatureSink(
            self.OUT_FC,
            'Output Strike/Dip Layer',
            type=QgsProcessing.TypeVectorPoint,
        ))

        # 2 – Elevation raster
        self.addParameter(QgsProcessingParameterRasterLayer(
            self.IN_RASTER,
            'Elevation Raster',
        ))

        # 3 – Triplet ID field (optional)
        p = QgsProcessingParameterField(
            self.TRIPLET_FIELD,
            'Triplet ID Field (optional — groups by feature order if omitted)',
            parentLayerParameterName=self.IN_POINTS,
            optional=True,
        )
        self.addParameter(p)

    # -------------------------------------------------------------------------
    # Main execution — mirrors execute() in the arcpy version
    # -------------------------------------------------------------------------

    def processAlgorithm(self, parameters, context, feedback):
        source        = self.parameterAsSource(parameters, self.IN_POINTS, context)
        raster_layer  = self.parameterAsRasterLayer(parameters, self.IN_RASTER, context)
        triplet_field = self.parameterAsString(parameters, self.TRIPLET_FIELD, context) or None

        if source is None:
            raise QgsProcessingException('Invalid input point layer.')
        if QgsWkbTypes.geometryType(source.wkbType()) != QgsWkbTypes.PointGeometry:
            raise QgsProcessingException('The Point Triplets Layer must be a point feature layer.')

        # Build output fields (mirrors out_fields in the arcpy version)
        out_fields = QgsFields()
        out_fields.append(QgsField('STRIKE',   QVariant.Double))
        out_fields.append(QgsField('DIP',      QVariant.Double))
        out_fields.append(QgsField('SRC_OID1', QVariant.LongLong))
        out_fields.append(QgsField('SRC_OID2', QVariant.LongLong))
        out_fields.append(QgsField('SRC_OID3', QVariant.LongLong))

        sink, dest_id = self.parameterAsSink(
            parameters, self.OUT_FC, context,
            out_fields, QgsWkbTypes.Point, source.sourceCrs(),
        )
        if sink is None:
            raise QgsProcessingException('Could not create output layer.')

        src_crs = source.sourceCrs()
        transform_context = context.transformContext()

        # Open DEM for elevation sampling (replaces _sample_elevations)
        dem_path = raster_layer.source()
        dem_ds = gdal.Open(dem_path)
        if dem_ds is None:
            raise QgsProcessingException(f'Cannot open elevation raster: {dem_path}')

        # Group features into triplets
        feedback.pushInfo('Reading input features...')
        if triplet_field:
            triplet_map = {}
            for feat in source.getFeatures():
                tid = feat[triplet_field]
                if tid is None:
                    feedback.pushWarning(f'Feature FID {feat.id()} has null TripletID; skipping.')
                    continue
                triplet_map.setdefault(tid, []).append(feat)
            bad = [tid for tid, feats in triplet_map.items() if len(feats) != 3]
            if bad:
                raise QgsProcessingException(
                    f'TripletID groups without exactly 3 points: '
                    f'{bad[:10]}{"..." if len(bad) > 10 else ""}'
                )
            groups = [sorted(feats, key=lambda f: f.id()) for feats in triplet_map.values()]
        else:
            all_feats = sorted(source.getFeatures(), key=lambda f: f.id())
            if len(all_feats) % 3 != 0:
                raise QgsProcessingException(
                    f'Input has {len(all_feats)} points; not divisible by 3. '
                    'Add/remove points or specify a Triplet ID field.'
                )
            groups = list(_batched(all_feats, 3))

        feedback.pushInfo(f'Found {len(groups)} triplets in input.')

        inserted = 0
        skipped = 0

        for i, trip in enumerate(groups):
            if feedback.isCanceled():
                break
            feedback.setProgress(int(100 * i / len(groups)))

            try:
                # Sample elevation from DEM for each of the three points
                items = []     # [(fid, x, y, z), ...]
                for feat in trip:
                    pt = feat.geometry().asPoint()
                    z = _sample_raster_at_xy(dem_ds, pt.x(), pt.y())
                    if z is None:
                        raise RuntimeError(
                            f'Could not sample DEM at FID {feat.id()} '
                            f'({pt.x():.2f}, {pt.y():.2f}).'
                        )
                    items.append((feat.id(), pt.x(), pt.y(), z))

                # Sort by FID ascending (mirrors items_sorted in arcpy version)
                items_sorted = sorted(items, key=lambda v: v[0])

                # Output point at location of middle-OID point (index 1)
                cx = items_sorted[1][1]
                cy = items_sorted[1][2]

                # Compute strike and dip from 3D points
                p1 = items_sorted[0][1:4]
                p2 = items_sorted[1][1:4]
                p3 = items_sorted[2][1:4]
                strike, dip = _strike_dip_from_three_points(p1, p2, p3)
                if math.isnan(strike) or math.isnan(dip):
                    feedback.pushWarning(
                        f'Degenerate triplet (collinear/coincident); '
                        f'skipping FIDs {[v[0] for v in items_sorted]}.'
                    )
                    continue

                out_feat = QgsFeature(out_fields)
                out_feat.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(cx, cy)))
                out_feat['STRIKE']   = strike
                out_feat['DIP']      = dip
                out_feat['SRC_OID1'] = items_sorted[0][0]
                out_feat['SRC_OID2'] = items_sorted[1][0]
                out_feat['SRC_OID3'] = items_sorted[2][0]
  
                sink.addFeature(out_feat, QgsFeatureSink.FastInsert)
                inserted += 1

            except Exception as ex:
                feedback.pushWarning(f'Failed triplet at index {i}: {ex}')

        dem_ds = None
        feedback.setProgress(100)
        feedback.pushInfo(
            f'Inserted {inserted} new strike/dip point(s). '
            f'Skipped {skipped} duplicate(s) within this run.'
        )
        return {self.OUT_FC: dest_id}


# =============================================================================
# Algorithm 2: Project Plane Across DEM (Outcrop Trace)
# =============================================================================

class ProjectPlaneAcrossDEMAlgorithm(QgsProcessingAlgorithm):
    """
    For each input point with a known strike/dip, computes the theoretical
    outcrop trace — the line where the geological plane intersects the DEM
    surface — within a specified radius.

    Direct port of the ProjectPlaneAcrossDEM tool from ThreePointProblemToolbox.pyt.
    arcpy.sa.Contour is replaced by gdal.ContourGenerate.
    arcpy.sa.ExtractByMask is replaced by gdal.Warp with a cutline shapefile.
    """

    IN_POINTS    = 'in_points'
    USE_SELECTED = 'use_selected'
    STRIKE_FIELD = 'strike_field'
    DIP_FIELD    = 'dip_field'
    IN_DEM       = 'in_dem'
    RESOLUTION   = 'resolution'
    RADIUS       = 'radius'
    OUT_FC       = 'out_fc'

    # -------------------------------------------------------------------------
    # Algorithm identity
    # -------------------------------------------------------------------------

    def createInstance(self):
        return ProjectPlaneAcrossDEMAlgorithm()

    def name(self):
        return 'projectplaneacrossdem'

    def displayName(self):
        return 'Project Plane Across DEM (Outcrop Trace)'

    def group(self):
        return 'Structural Geology'

    def groupId(self):
        return 'structuralgeology'

    def shortHelpString(self):
        return (
            'For each input point with a known strike/dip, computes the theoretical '
            'outcrop trace — the line where the geological plane intersects the DEM '
            'surface — within a specified radius.\n\n'
            'For each point the tool:\n'
            '  1. Clips the DEM to a circular buffer of the given radius.\n'
            '  2. Builds a plane elevation grid anchored to the sampled Z at the '
            'point centre: z_plane = z0 − tan(dip) × (distance along dip direction).\n'
            '  3. Computes a difference raster: diff = DEM − plane.\n'
            '  4. Extracts the zero-elevation contour (where DEM = plane) using '
            'gdal.ContourGenerate as the outcrop trace.\n'
            '  5. Clips the trace to the circle and appends polyline segments to the '
            'output layer.\n\n'
            'Strike/dip fields are auto-detected from common names if not specified.\n\n'
            'Note: "Use selected features only" mirrors the arcpy version. In QGIS '
            'you can also control this via the Processing dialog\'s feature filter.'
        )

    # -------------------------------------------------------------------------
    # Parameter definitions — mirrors getParameterInfo() in the arcpy version
    # -------------------------------------------------------------------------

    def initAlgorithm(self, config=None):
        # 0 – Input points with strike/dip
        self.addParameter(QgsProcessingParameterFeatureSource(
            self.IN_POINTS,
            'Input Points (with Strike/Dip)',
            [QgsProcessing.TypeVectorPoint],
        ))

        # 1 – Use selected features only
        self.addParameter(QgsProcessingParameterBoolean(
            self.USE_SELECTED,
            'Use selected features only',
            defaultValue=False,
            optional=True,
        ))

        # 2 – Strike field (optional, auto-detected)
        self.addParameter(QgsProcessingParameterField(
            self.STRIKE_FIELD,
            'Strike Field (optional — auto-detected if omitted)',
            parentLayerParameterName=self.IN_POINTS,
            type=QgsProcessingParameterField.Numeric,
            optional=True,
        ))

        # 3 – Dip field (optional, auto-detected)
        self.addParameter(QgsProcessingParameterField(
            self.DIP_FIELD,
            'Dip Field (optional — auto-detected if omitted)',
            parentLayerParameterName=self.IN_POINTS,
            type=QgsProcessingParameterField.Numeric,
            optional=True,
        ))

        # 4 – DEM raster
        self.addParameter(QgsProcessingParameterRasterLayer(
            self.IN_DEM,
            'DEM Raster',
        ))

        # 5 – Target resolution (optional)
        p = QgsProcessingParameterNumber(
            self.RESOLUTION,
            'Target Resolution (map units; optional — uses DEM native resolution if omitted)',
            type=QgsProcessingParameterNumber.Double,
            optional=True,
            minValue=0.0,
        )
        self.addParameter(p)

        # 6 – Radius
        self.addParameter(QgsProcessingParameterNumber(
            self.RADIUS,
            'Radius (map units)',
            type=QgsProcessingParameterNumber.Double,
            minValue=0.0,
        ))

        # 7 – Output polyline feature class
        self.addParameter(QgsProcessingParameterFeatureSink(
            self.OUT_FC,
            'Output Outcrop Trace Polylines',
            type=QgsProcessing.TypeVectorLine,
        ))

    # -------------------------------------------------------------------------
    # Main execution — mirrors execute() in the arcpy version
    # -------------------------------------------------------------------------

    def processAlgorithm(self, parameters, context, feedback):
        source       = self.parameterAsSource(parameters, self.IN_POINTS, context)
        use_selected = self.parameterAsBoolean(parameters, self.USE_SELECTED, context)
        strike_field = self.parameterAsString(parameters, self.STRIKE_FIELD, context) or None
        dip_field    = self.parameterAsString(parameters, self.DIP_FIELD, context) or None
        dem_layer    = self.parameterAsRasterLayer(parameters, self.IN_DEM, context)
        resolution   = (self.parameterAsDouble(parameters, self.RESOLUTION, context)
                        if parameters.get(self.RESOLUTION) is not None else None)
        radius       = self.parameterAsDouble(parameters, self.RADIUS, context)

        if source is None:
            raise QgsProcessingException('Invalid input point layer.')
        if QgsWkbTypes.geometryType(source.wkbType()) != QgsWkbTypes.PointGeometry:
            raise QgsProcessingException('Input points must be a point feature layer.')

        # Auto-detect strike/dip fields (mirrors _autodetect_fields in arcpy version)
        strike_field = _autodetect_field(
            source, strike_field, ['strike', 'strike_deg', 'strk']
        )
        dip_field = _autodetect_field(
            source, dip_field, ['dip', 'dip_deg']
        )
        if not strike_field:
            raise QgsProcessingException(
                'Strike field not provided and could not be auto-detected. '
                'Expected: strike, strike_deg, or strk.'
            )
        if not dip_field:
            raise QgsProcessingException(
                'Dip field not provided and could not be auto-detected. '
                'Expected: dip or dip_deg.'
            )

        # DEM metadata
        dem_path = dem_layer.source()
        dem_info_ds = gdal.Open(dem_path)
        if dem_info_ds is None:
            raise QgsProcessingException(f'Cannot open DEM: {dem_path}')
        dem_gt   = dem_info_ds.GetGeoTransform()
        dem_proj = dem_info_ds.GetProjection()
        base_cs  = dem_gt[1]   # native cell width
        dem_info_ds = None

        target_cs = resolution if (resolution and resolution > 0) else base_cs
        feedback.pushInfo(f'[Setup] radius={radius}, target_cs={target_cs}')

        # Build output fields (SRC_OID, STRIKE, DIP — matches arcpy version)
        out_fields = QgsFields()
        out_fields.append(QgsField('SRC_OID', QVariant.LongLong))
        out_fields.append(QgsField('STRIKE',  QVariant.Double))
        out_fields.append(QgsField('DIP',     QVariant.Double))

        sink, dest_id = self.parameterAsSink(
            parameters, self.OUT_FC, context,
            out_fields, QgsWkbTypes.LineString, dem_layer.crs(),
        )
        if sink is None:
            raise QgsProcessingException('Could not create output layer.')

        # Collect features, optionally filtering to selection
        all_feats = list(source.getFeatures())
        if use_selected:
            # Attempt to access the map layer for its selection set
            layer_id = parameters.get(self.IN_POINTS)
            map_layer = context.getMapLayer(layer_id) if isinstance(layer_id, str) else None
            if map_layer is not None:
                sel_ids = set(map_layer.selectedFeatureIds())
                if sel_ids:
                    all_feats = [f for f in all_feats if f.id() in sel_ids]
                    if not all_feats:
                        feedback.pushWarning(
                            'use_selected=True but no features are currently selected; '
                            'processing all features.'
                        )
                        all_feats = list(source.getFeatures())
            else:
                feedback.pushWarning(
                    'use_selected=True but could not access map layer; '
                    'processing all features.'
                )

        npts = len(all_feats)
        feedback.pushInfo(f'[Setup] {npts} points to process.')

        # CRS reprojection: points → DEM CRS if they differ
        src_crs    = source.sourceCrs()
        dem_crs    = dem_layer.crs()
        xform      = None
        if src_crs.authid() != dem_crs.authid():
            xform = QgsCoordinateTransform(src_crs, dem_crs, context.transformContext())
            feedback.pushInfo('[Setup] Projecting points to DEM CRS.')

        inserted_total = 0
        tmpdir = tempfile.mkdtemp(prefix='qgis_3pt_')

        try:
            for i, feat in enumerate(all_feats):
                if feedback.isCanceled():
                    break
                feedback.setProgress(int(100 * i / max(npts, 1)))

                src_oid = feat.id()
                pt = feat.geometry().asPoint()

                try:
                    strike = float(feat[strike_field])
                    dip    = float(feat[dip_field])
                except Exception:
                    feedback.pushWarning(
                        f'[Point FID {src_oid}] Strike/Dip conversion failed; skipping.'
                    )
                    continue

                # Reproject to DEM CRS if needed
                px, py = pt.x(), pt.y()
                if xform is not None:
                    p_dem = xform.transform(QgsPointXY(px, py))
                    px, py = p_dem.x(), p_dem.y()

                feedback.pushInfo(
                    f'\n[Point {i + 1}/{npts}] FID={src_oid} strike={strike} dip={dip}'
                )
                t0 = time.time()

                try:
                    # ── Build circle shapefile for GDAL cutline ──────────────
                    circle_shp = os.path.join(tmpdir, f'circle_{src_oid}.shp')
                    circle_geom = _write_circle_shapefile(px, py, radius, dem_proj, circle_shp)

                    # ── Clip DEM to circle (replaces arcpy.sa.ExtractByMask) ─
                    clipped_tif = os.path.join(tmpdir, f'dem_clip_{src_oid}.tif')
                    warp_opts = gdal.WarpOptions(
                        cutlineDSName=circle_shp,
                        cropToCutline=True,
                        dstNodata=-9999.0,
                        xRes=target_cs,
                        yRes=target_cs,
                        resampleAlg=gdal.GRA_Bilinear,
                        outputType=gdal.GDT_Float32,
                    )
                    gdal.Warp(clipped_tif, dem_path, options=warp_opts)

                    clip_ds = gdal.Open(clipped_tif)
                    if clip_ds is None:
                        feedback.pushWarning(
                            f'[Point FID {src_oid}] DEM clip failed; skipping.'
                        )
                        continue

                    clip_gt  = clip_ds.GetGeoTransform()
                    ncols    = clip_ds.RasterXSize
                    nrows    = clip_ds.RasterYSize
                    csx      = clip_gt[1]
                    csy      = abs(clip_gt[5])
                    x_left   = clip_gt[0]
                    y_top    = clip_gt[3]

                    feedback.pushInfo(
                        f'[Point FID {src_oid}] clip: '
                        f'({x_left:.1f}, {y_top - nrows * csy:.1f}) – '
                        f'({x_left + ncols * csx:.1f}, {y_top:.1f}), '
                        f'cols={ncols}, rows={nrows}, cs={csx:.3f}'
                    )

                    # ── Sample z0 at point centre ────────────────────────────
                    z0 = _sample_raster_at_xy(clip_ds, px, py)
                    if z0 is None:
                        feedback.pushWarning(
                            f'[Point FID {src_oid}] Could not sample DEM at centre; skipping.'
                        )
                        clip_ds = None
                        continue
                    feedback.pushInfo(f'[Point FID {src_oid}] z0={z0:.3f}')

                    # ── Build plane elevation grid (portable numpy math) ─────
                    cols = np.arange(ncols, dtype=float)
                    rows = np.arange(nrows, dtype=float)
                    Xc = x_left + (cols + 0.5) * csx
                    Yc = y_top  - (rows + 0.5) * csy
                    Xg, Yg = np.meshgrid(Xc, Yc)      # shape (nrows, ncols)

                    dipdir = (strike + 90.0) % 360.0
                    diprad = math.radians(dip)
                    dd_rad = math.radians(dipdir)
                    ux     = math.sin(dd_rad)
                    uy     = math.cos(dd_rad)
                    slope  = math.tan(diprad)

                    # Signed distance along dip direction from the centre
                    s       = (Xg - px) * ux + (Yg - py) * uy
                    z_plane = z0 - slope * s

                    # ── Difference raster: DEM − plane ───────────────────────
                    clip_band = clip_ds.GetRasterBand(1)
                    dem_arr   = clip_band.ReadAsArray().astype(float)
                    nodata    = clip_band.GetNoDataValue()
                    if nodata is not None:
                        dem_arr[dem_arr == nodata] = np.nan

                    diff_arr = dem_arr - z_plane
                    diff_arr = np.where(np.isnan(diff_arr), -9999.0, diff_arr)
                    clip_ds  = None   # close clipped DEM

                    # Write diff to temp GeoTIFF
                    diff_tif  = os.path.join(tmpdir, f'diff_{src_oid}.tif')
                    drv       = gdal.GetDriverByName('GTiff')
                    diff_ds   = drv.Create(diff_tif, ncols, nrows, 1, gdal.GDT_Float32)
                    diff_ds.SetGeoTransform(clip_gt)
                    diff_ds.SetProjection(dem_proj)
                    diff_band = diff_ds.GetRasterBand(1)
                    diff_band.WriteArray(diff_arr.astype(np.float32))
                    diff_band.SetNoDataValue(-9999.0)
                    diff_ds.FlushCache()

                    # ── Extract 0-contour (replaces arcpy.sa.Contour) ────────
                    # interval=1000000, base=0 → only the 0-level contour is
                    # ever present in a local-area diff raster.
                    mem_driver   = ogr.GetDriverByName('Memory')
                    contour_ds   = mem_driver.CreateDataSource('')
                    contour_srs  = osr.SpatialReference()
                    contour_srs.ImportFromWkt(dem_proj)
                    contour_lyr  = contour_ds.CreateLayer(
                        'contours', srs=contour_srs, geom_type=ogr.wkbLineString
                    )
                    contour_lyr.CreateField(ogr.FieldDefn('ID',   ogr.OFTInteger))
                    contour_lyr.CreateField(ogr.FieldDefn('ELEV', ogr.OFTReal))

                    gdal.ContourGenerate(
                        diff_band,      # source band
                        1000000.0,      # contour interval
                        0.0,            # base contour
                        [],             # fixed levels (none; use interval)
                        1,              # useNoData
                        -9999.0,        # noDataValue
                        contour_lyr,    # destination OGR layer
                        0,              # idField index
                        1,              # elevField index
                    )
                    diff_ds = None   # close diff raster

                    # ── Clip contours to circle and write to sink ────────────
                    n_added = 0
                    contour_lyr.ResetReading()
                    for c_feat in contour_lyr:
                        c_geom = c_feat.GetGeometryRef()
                        if c_geom is None:
                            continue
                        clipped_geom = c_geom.Intersection(circle_geom)
                        if clipped_geom is None or clipped_geom.IsEmpty():
                            continue
                        qgs_geom = QgsGeometry.fromWkt(clipped_geom.ExportToWkt())
                        if qgs_geom.isEmpty():
                            continue
                        out_feat = QgsFeature(out_fields)
                        out_feat.setGeometry(qgs_geom)
                        out_feat['SRC_OID'] = int(src_oid)
                        out_feat['STRIKE']  = float(strike)
                        out_feat['DIP']     = float(dip)
                        sink.addFeature(out_feat, QgsFeatureSink.FastInsert)
                        n_added += 1

                    contour_ds = None
                    inserted_total += n_added

                    if n_added == 0:
                        feedback.pushWarning(
                            f'[Point FID {src_oid}] No contour segments produced. '
                            'Try increasing radius or check DEM coverage.'
                        )
                    else:
                        feedback.pushInfo(
                            f'[Point FID {src_oid}] appended {n_added} segment(s). '
                            f'({time.time() - t0:.1f}s)'
                        )

                except Exception as ex:
                    feedback.pushWarning(f'[Point FID {src_oid}] Failed: {ex}')
                finally:
                    gc.collect()

        finally:
            # Clean up all temp files written during this run
            import shutil
            try:
                shutil.rmtree(tmpdir, ignore_errors=True)
            except Exception:
                pass

        feedback.setProgress(100)
        feedback.pushInfo(
            f'[Finish] Inserted {inserted_total} segment(s) across {npts} point(s).'
        )
        return {self.OUT_FC: dest_id}

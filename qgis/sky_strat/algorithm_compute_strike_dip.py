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

from .shared import _sample_raster_at_xy

gdal.UseExceptions()

def _strike_dip_from_three_points(p1, p2, p3):
    """
    p1, p2, p3: (x, y, z). Returns (strike_deg, dip_deg) using RHR convention.
    Normal forced to point upward (nz > 0). Returns (nan, nan) for degenerate inputs.
    """
    v1 = (p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2])
    v2 = (p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2])

    nx = v1[1] * v2[2] - v1[2] * v2[1]
    ny = v1[2] * v2[0] - v1[0] * v2[2]
    nz = v1[0] * v2[1] - v1[1] * v2[0]

    # Check for degenerate case of collinear or coincident points
    if math.sqrt(nx*nx + ny*ny + nz*nz) == 0:
        return (float("nan"), float("nan"))

    if nz < 0:
        nx, ny, nz = -nx, -ny, -nz

    h = math.sqrt(nx*nx + ny*ny)
    dip_deg = math.degrees(math.atan2(h, nz))
    dipdir = math.degrees(math.atan2(nx, ny))
    if dipdir < 0:
        dipdir += 360.0
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

    # ------------------------
    # Parameter definitions
    # ------------------------

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
    # Main execution
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
                items = [] # [(fid, x, y, z), ...]
                for j, feat in enumerate(trip):
                    pt = feat.geometry().asPoint()
                    z = _sample_raster_at_xy(dem_ds, pt.x(), pt.y())
                    if z is None:
                        raise RuntimeError(
                            f'Could not sample DEM at FID {feat.id()} '
                            f'({pt.x():.2f}, {pt.y():.2f}).'
                        )
                    feedback.pushInfo(f'XYZ for point {j} in triplet {i}: ({pt.x():.2f}, {pt.y():.2f}, {z:.2f})')
                    items.append((feat.id(), pt.x(), pt.y(), z))

                # Sort by FID ascending (mirrors items_sorted in arcpy version)
                items_sorted = sorted(items, key=lambda v: v[0])

                # Output point at location of closest to centroid
                centroid_x = (items_sorted[1][1] + items_sorted[2][1] + items_sorted[0][1])/3.0
                centroid_y = (items_sorted[1][2] + items_sorted[2][2] + items_sorted[0][2])/3.0

                nearest_idx = min(range(3), key=lambda j: math.hypot(items_sorted[j][1] - centroid_x, items_sorted[j][2] - centroid_y))
                cx, cy = items_sorted[nearest_idx][1], items_sorted[nearest_idx][2]
                # Compute strike and dip from 3D points
                p1 = items_sorted[0][1:4] # x,y,z
                p2 = items_sorted[1][1:4]
                p3 = items_sorted[2][1:4]
                strike, dip = _strike_dip_from_three_points(p1, p2, p3)
                if math.isnan(strike) or math.isnan(dip):
                    feedback.pushWarning(
                        f'Degenerate triplet (collinear/coincident); '
                        f'skipping FIDs {[v[0] for v in items_sorted]}.'
                    )
                    continue
                feedback.pushInfo(f'Strike/dip {j} in triplet {i}: {strike:.2f}/{dip:.2f}')

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
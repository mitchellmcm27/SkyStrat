import math
import os
import tempfile

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from osgeo import gdal, ogr, osr

from qgis.core import (
    QgsFeature,
    QgsFeatureSink,
    QgsField,
    QgsFields,
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingException,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterFeatureSink,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterField,
    QgsProcessingParameterFileDestination,
    QgsProcessingParameterNumber,
    QgsProcessingParameterRasterDestination,
    QgsProcessingParameterRasterLayer,
    QgsWkbTypes,
)
from qgis.PyQt.QtCore import (QCoreApplication, QVariant)

gdal.UseExceptions()


class BuskDownPlungeAlgorithm(QgsProcessingAlgorithm):
    """QGIS Processing Algorithm implementing the Busk downplunge view method."""

    # Parameter name constants (mirror arcpy parameter names)
    STRIKE_DIP_FC              = 'strike_dip_fc'
    STRIKE_FIELD               = 'strike_field'
    DIP_FIELD                  = 'dip_field'
    OVERTURNED_FIELD           = 'overturned_field'
    LABEL_FIELD                = 'label_field'
    STRAT_HEIGHT_FIELD         = 'strat_height_field'
    DEM_RASTER                 = 'dem_raster'
    STUDY_AREA                 = 'study_area'
    FOLD_TREND                 = 'fold_trend'
    FOLD_PLUNGE                = 'fold_plunge'
    OUTPUT_PDF                 = 'output_pdf'
    PLOT_DEM_CELLS             = 'plot_dem_cells'
    WEDGE_RASTER_OUTPUT        = 'wedge_raster_output'
    OUTPUT_STRIKE_DIP_FC       = 'output_strike_dip_fc'
    STRAT_HEIGHT_RASTER_OUTPUT = 'strat_height_raster_output'

    # -------------------------------------------------------------------------
    # Algorithm identity
    # -------------------------------------------------------------------------

    def createInstance(self):
        return BuskDownPlungeAlgorithm()

    def name(self):
        return 'buskdownplungeview'

    def displayName(self):
        return 'Stratigraphic Heights (Busk Method)'

    def shortHelpString(self):
        return (
            'Creates a downplunge projection of strike/dip measurements using the '
            'Busk Method (cylindrical fold geometry) and optionally calculates '
            'stratigraphic heights for measurement points and DEM cells.\n\n'
            'All inputs must use a projected coordinate system (UTM recommended).'
        )

    # -------------------------------------------------------------------------
    # Parameter definitions  (mirrors getParameterInfo in the arcpy version)
    # -------------------------------------------------------------------------

    def initAlgorithm(self, config=None):
        # 0 – Input strike/dip point feature class
        self.addParameter(QgsProcessingParameterFeatureSource(
            self.STRIKE_DIP_FC,
            'Strike and Dip Point Feature Class',
            [QgsProcessing.TypeVectorPoint],
        ))

        # 1 – Strike field
        self.addParameter(QgsProcessingParameterField(
            self.STRIKE_FIELD,
            'Strike Field',
            parentLayerParameterName=self.STRIKE_DIP_FC,
            type=QgsProcessingParameterField.Numeric,
        ))

        # 2 – Dip field
        self.addParameter(QgsProcessingParameterField(
            self.DIP_FIELD,
            'Dip Field',
            parentLayerParameterName=self.STRIKE_DIP_FC,
            type=QgsProcessingParameterField.Numeric,
        ))

        # 3 – Overturned flag field (optional)
        p = QgsProcessingParameterField(
            self.OVERTURNED_FIELD,
            'Overturned Flag Field (Optional — 1 = overturned)',
            parentLayerParameterName=self.STRIKE_DIP_FC,
            optional=True,
        )
        self.addParameter(p)

        # 4 – Label field
        p = QgsProcessingParameterField(
            self.LABEL_FIELD,
            'Label Field (for point identification)',
            parentLayerParameterName=self.STRIKE_DIP_FC,
        )
        p.setDefaultValue('OBJECTID')
        self.addParameter(p)

        # 5 – Stratigraphic height field (optional)
        p = QgsProcessingParameterField(
            self.STRAT_HEIGHT_FIELD,
            'Stratigraphic Height Field (Optional)',
            parentLayerParameterName=self.STRIKE_DIP_FC,
            type=QgsProcessingParameterField.Numeric,
            optional=True,
        )
        self.addParameter(p)

        # 6 – DEM raster
        self.addParameter(QgsProcessingParameterRasterLayer(
            self.DEM_RASTER,
            'Digital Elevation Model (DEM)',
        ))

        # 7 – Study area polygon
        self.addParameter(QgsProcessingParameterFeatureSource(
            self.STUDY_AREA,
            'Study Area Polygon',
            [QgsProcessing.TypeVectorPolygon],
        ))

        # 8 – Fold axis trend (optional)
        p = QgsProcessingParameterNumber(
            self.FOLD_TREND,
            'Fold Axis Trend (azimuth, optional — auto-calculated if omitted)',
            type=QgsProcessingParameterNumber.Double,
            optional=True,
            minValue=0.0,
            maxValue=360.0,
        )
        self.addParameter(p)

        # 9 – Fold axis plunge (optional)
        p = QgsProcessingParameterNumber(
            self.FOLD_PLUNGE,
            'Fold Axis Plunge (degrees, optional — auto-calculated if omitted)',
            type=QgsProcessingParameterNumber.Double,
            optional=True,
            minValue=0.0,
            maxValue=90.0,
        )
        self.addParameter(p)

        # 10 – Output PDF
        self.addParameter(QgsProcessingParameterFileDestination(
            self.OUTPUT_PDF,
            'Output PDF Path',
            fileFilter='PDF Files (*.pdf)',
        ))

        # 11 – Plot DEM cell centroids
        p = QgsProcessingParameterBoolean(
            self.PLOT_DEM_CELLS,
            'Plot DEM Cell Centroids on Profile',
            defaultValue=False,
            optional=True,
        )
        self.addParameter(p)

        # 12 – Wedge assignment raster (optional output)
        p = QgsProcessingParameterRasterDestination(
            self.WEDGE_RASTER_OUTPUT,
            'Output Wedge Assignment Raster (Optional)',
            optional=True,
            createByDefault=False,
        )
        self.addParameter(p)

        # 13 – Output feature class with heights (optional output)
        p = QgsProcessingParameterFeatureSink(
            self.OUTPUT_STRIKE_DIP_FC,
            'Output Strike/Dip Feature Class with Heights (Optional)',
            optional=True,
            createByDefault=False,
        )
        self.addParameter(p)

        # 14 – Stratigraphic height raster (optional output)
        p = QgsProcessingParameterRasterDestination(
            self.STRAT_HEIGHT_RASTER_OUTPUT,
            'Output Stratigraphic Height Raster (Optional)',
            optional=True,
            createByDefault=False,
        )
        self.addParameter(p)

    # -------------------------------------------------------------------------
    # Main execution  (mirrors execute() in the arcpy version)
    # -------------------------------------------------------------------------

    def processAlgorithm(self, parameters, context, feedback):
        # Read parameters
        strike_dip_source    = self.parameterAsSource(parameters, self.STRIKE_DIP_FC, context)
        strike_field         = self.parameterAsString(parameters, self.STRIKE_FIELD, context)
        dip_field            = self.parameterAsString(parameters, self.DIP_FIELD, context)
        overturned_field     = self.parameterAsString(parameters, self.OVERTURNED_FIELD, context) or None
        label_field          = self.parameterAsString(parameters, self.LABEL_FIELD, context)
        strat_height_field   = self.parameterAsString(parameters, self.STRAT_HEIGHT_FIELD, context) or None
        dem_layer            = self.parameterAsRasterLayer(parameters, self.DEM_RASTER, context)
        study_area_source    = self.parameterAsSource(parameters, self.STUDY_AREA, context)
        output_pdf           = self.parameterAsFileOutput(parameters, self.OUTPUT_PDF, context)
        plot_dem_cells       = self.parameterAsBoolean(parameters, self.PLOT_DEM_CELLS, context)

        fold_trend  = (self.parameterAsDouble(parameters, self.FOLD_TREND, context)
                       if parameters.get(self.FOLD_TREND) is not None else None)
        fold_plunge = (self.parameterAsDouble(parameters, self.FOLD_PLUNGE, context)
                       if parameters.get(self.FOLD_PLUNGE) is not None else None)

        wedge_raster_output        = parameters.get(self.WEDGE_RASTER_OUTPUT)
        output_strike_dip_fc_param = parameters.get(self.OUTPUT_STRIKE_DIP_FC)
        strat_height_raster_output = parameters.get(self.STRAT_HEIGHT_RASTER_OUTPUT)

        dem_path          = dem_layer.source()
        study_area_path   = self._save_source_to_shapefile(study_area_source, 'study_area', feedback)

        feedback.pushInfo('=' * 60)
        feedback.pushInfo('Busk Method - Downplunge View Generation')
        feedback.pushInfo('=' * 60)

        # Check coordinate systems
        self.check_utm_projection(strike_dip_source.sourceCrs(), 'Strike/Dip Feature Class', feedback)
        self.check_utm_projection(dem_layer.crs(), 'DEM Raster', feedback)
        self.check_utm_projection(study_area_source.sourceCrs(), 'Study Area Polygon', feedback)

        # Read strike and dip data
        feedback.pushInfo('\nReading strike and dip measurements...')
        strike_dip_data = self.read_strike_dip_data(
            strike_dip_source, strike_field, dip_field,
            overturned_field, label_field, strat_height_field,
            dem_path, feedback,
        )
        feedback.pushInfo(f"  Loaded {len(strike_dip_data['x'])} measurements")

        # Calculate or use provided fold axis
        if fold_trend is None or fold_plunge is None:
            feedback.pushInfo('\nCalculating best-fit fold axis...')
            fold_trend, fold_plunge = self.calculate_fold_axis(strike_dip_data)
            feedback.pushInfo(f'  Calculated fold axis: {fold_trend:.1f}° / {fold_plunge:.1f}°')
        else:
            feedback.pushInfo(f'\nUsing provided fold axis: {fold_trend:.1f}° / {fold_plunge:.1f}°')

        # Calculate projected attitudes on profile plane
        feedback.pushInfo('\nCalculating projected attitudes on profile plane...')
        strike_dip_data = self.calculate_projected_attitudes(strike_dip_data, fold_trend, fold_plunge)
        feedback.pushInfo(f"  Projected {len(strike_dip_data['x'])} attitudes")

        # Sort points by profile x-coordinate
        profile_x_coords = strike_dip_data['profile_x']
        sorted_indices   = np.argsort(profile_x_coords)

        # Log projected dip angles
        feedback.pushInfo('\nProjected dip angles in profile plane:')
        for i in range(len(strike_dip_data['labels'])):
            att        = strike_dip_data['projected_attitudes'][i]
            dip_angle  = np.degrees(np.arctan2(abs(att['y']), abs(att['x'])))
            apparent   = 90 - dip_angle if abs(att['y']) > abs(att['x']) else dip_angle
            feedback.pushInfo(f"  {strike_dip_data['labels'][i]}: {apparent:.1f}°")

        # Log stratigraphic height vectors
        feedback.pushInfo('\nStratigraphic height vectors in profile plane (x=right, y=up):')
        for i in range(len(strike_dip_data['labels'])):
            sv     = strike_dip_data['strat_height_vectors'][i]
            status = 'OVERTURNED' if strike_dip_data['overturned'][i] else 'UPRIGHT'
            feedback.pushInfo(
                f"  {strike_dip_data['labels'][i]} ({status}): "
                f"x={sv['x']:+.4f}, y={sv['y']:+.4f}"
            )

        # Analyze wedges
        feedback.pushInfo('\nAnalyzing wedge geometry between adjacent points...')
        wedge_data = self.analyze_wedges(strike_dip_data, sorted_indices)

        # Log parallel bed info
        feedback.pushInfo('\nParallel bed analysis:')
        has_parallel = False
        for wedge in wedge_data:
            if wedge['type'] in ('parallel', 'invalid_parallel'):
                has_parallel = True
                idx1, idx2 = wedge['point1_idx'], wedge['point2_idx']
                sv1 = strike_dip_data['strat_height_vectors'][idx1]
                sv2 = strike_dip_data['strat_height_vectors'][idx2]
                feedback.pushInfo(
                    f"  {strike_dip_data['labels'][idx1]} and "
                    f"{strike_dip_data['labels'][idx2]}: PARALLEL beds"
                )
                feedback.pushInfo(f"    {strike_dip_data['labels'][idx1]} strat height vector: x={sv1['x']:+.4f}, y={sv1['y']:+.4f}")
                feedback.pushInfo(f"    {strike_dip_data['labels'][idx2]} strat height vector: x={sv2['x']:+.4f}, y={sv2['y']:+.4f}")
                if wedge['type'] == 'parallel':
                    feedback.pushInfo('    -> VALID: Consistent stratigraphic directions')
                else:
                    feedback.pushInfo('    -> INVALID: Opposite stratigraphic directions')
        if not has_parallel:
            feedback.pushInfo('  No parallel beds found')
        feedback.pushInfo(f'  Analyzed {len(wedge_data)} wedges')

        # Extract DEM cell centroids if needed
        dem_cell_data = None
        if plot_dem_cells or wedge_raster_output:
            if not plot_dem_cells:
                feedback.pushInfo('  Extracting DEM cell centroids for wedge assignment raster...')
            else:
                feedback.pushInfo('  Extracting DEM cell centroids...')
            dem_cell_data = self.extract_dem_cell_centroids(dem_path, study_area_path, feedback)
            feedback.pushInfo(f"  Extracted {len(dem_cell_data['x'])} DEM cells")

        # Build output sink for feature class (if requested)
        output_sink      = None
        output_sink_id   = None
        output_fc_fields = None
        if output_strike_dip_fc_param not in (None, ''):
            # Build the field list now (need it before create_downplunge_view)
            output_fc_fields = QgsFields()
            for f in strike_dip_source.fields():
                output_fc_fields.append(f)
            output_fc_fields.append(QgsField('strat_height_out', QVariant.Double))
            output_sink, output_sink_id = self.parameterAsSink(
                parameters, self.OUTPUT_STRIKE_DIP_FC, context,
                output_fc_fields,
                strike_dip_source.wkbType(),
                strike_dip_source.sourceCrs(),
            )

        # Resolve output raster paths (empty string = not requested)
        wedge_raster_path = (
            self.parameterAsOutputLayer(parameters, self.WEDGE_RASTER_OUTPUT, context)
            if wedge_raster_output not in (None, '') else None
        )
        strat_raster_path = (
            self.parameterAsOutputLayer(parameters, self.STRAT_HEIGHT_RASTER_OUTPUT, context)
            if strat_height_raster_output not in (None, '') else None
        )

        # Create downplunge projection (main analysis + plotting)
        feedback.pushInfo('\nCreating downplunge projection...')
        self.create_downplunge_view(
            strike_dip_data, wedge_data, dem_cell_data,
            fold_trend, fold_plunge,
            output_pdf, dem_path, study_area_path,
            wedge_raster_path,
            strike_dip_source, output_sink, label_field, plot_dem_cells,
            strat_raster_path,
            feedback,
        )
        feedback.pushInfo(f'  Downplunge view saved to: {output_pdf}')

        feedback.pushInfo('\n' + '=' * 60)
        feedback.pushInfo('Processing complete!')
        feedback.pushInfo('=' * 60)

        results = {self.OUTPUT_PDF: output_pdf}
        if output_sink_id:
            results[self.OUTPUT_STRIKE_DIP_FC] = output_sink_id
        if wedge_raster_path:
            results[self.WEDGE_RASTER_OUTPUT] = wedge_raster_path
        if strat_raster_path:
            results[self.STRAT_HEIGHT_RASTER_OUTPUT] = strat_raster_path
        return results

    # =========================================================================
    # Utility helpers  (no arcpy equivalent in the original)
    # =========================================================================

    def _save_source_to_shapefile(self, source, name, feedback):
        """Write a QgsProcessingFeatureSource to a temporary shapefile for GDAL."""
        temp_dir  = tempfile.mkdtemp()
        temp_path = os.path.join(temp_dir, f'{name}.shp')

        srs = osr.SpatialReference()
        srs.ImportFromWkt(source.sourceCrs().toWkt())

        geom_type_map = {
            QgsWkbTypes.Point:        ogr.wkbPoint,
            QgsWkbTypes.PointZ:       ogr.wkbPoint25D,
            QgsWkbTypes.Polygon:      ogr.wkbPolygon,
            QgsWkbTypes.PolygonZ:     ogr.wkbPolygon25D,
            QgsWkbTypes.MultiPolygon: ogr.wkbMultiPolygon,
        }
        ogr_geom = geom_type_map.get(source.wkbType(), ogr.wkbUnknown)

        type_map = {
            QVariant.Int:      ogr.OFTInteger,
            QVariant.LongLong: ogr.OFTInteger64,
            QVariant.Double:   ogr.OFTReal,
            QVariant.String:   ogr.OFTString,
        }

        driver = ogr.GetDriverByName('ESRI Shapefile')
        ds     = driver.CreateDataSource(temp_path)
        lyr    = ds.CreateLayer(name, srs, ogr_geom)

        field_names = []
        for qf in source.fields():
            fname = qf.name()[:10]  # shapefile limit
            ogr_type = type_map.get(qf.type(), ogr.OFTString)
            lyr.CreateField(ogr.FieldDefn(fname, ogr_type))
            field_names.append((qf.name(), fname))

        for feat in source.getFeatures():
            ogr_feat = ogr.Feature(lyr.GetLayerDefn())
            ogr_feat.SetGeometry(ogr.CreateGeometryFromWkt(feat.geometry().asWkt()))
            for qname, fname in field_names:
                val = feat[qname]
                if val is not None:
                    try:
                        ogr_feat.SetField(fname, val)
                    except Exception:
                        pass
            lyr.CreateFeature(ogr_feat)

        ds = None
        return temp_path

    def _sample_raster(self, band, gt, nodata, x, y):
        """Sample a GDAL raster band at geographic coordinate (x, y)."""
        if band is None or gt is None:
            return None
        px = int((x - gt[0]) / gt[1])
        py = int((y - gt[3]) / gt[5])
        if px < 0 or py < 0 or px >= band.XSize or py >= band.YSize:
            return None
        val = band.ReadAsArray(px, py, 1, 1)
        if val is None:
            return None
        val = float(val[0][0])
        if nodata is not None and val == nodata:
            return None
        return val

    def _open_clipped_raster(self, dem_path, study_area_path):
        """Clip DEM to study area and return (gdal.Dataset, temp_path)."""
        temp_dir     = tempfile.mkdtemp()
        clipped_path = os.path.join(temp_dir, 'dem_clipped.tif')
        gdal.Warp(
            clipped_path, dem_path,
            cutlineDSName=study_area_path,
            cropToCutline=True,
            dstNodata=-9999,
            format='GTiff',
        )
        ds = gdal.Open(clipped_path)
        return ds, clipped_path

    # =========================================================================
    # check_utm_projection  (replaces arcpy.Describe + AddWarning/ValueError)
    # =========================================================================

    def check_utm_projection(self, crs, dataset_name, feedback):
        """Raise an error if CRS is geographic; warn if not UTM."""
        if crs.isGeographic():
            raise QgsProcessingException(
                f'{dataset_name} must be in a projected coordinate system (UTM required)'
            )
        if 'utm' not in crs.description().lower():
            feedback.pushWarning(
                f'WARNING: {dataset_name} may not be in UTM (found: {crs.description()})'
            )
            feedback.pushWarning('  The tool expects UTM coordinates for accurate calculations')
        else:
            feedback.pushInfo(f'  {dataset_name}: {crs.description()} ✓')

    # =========================================================================
    # read_strike_dip_data  (replaces arcpy.Describe + arcpy.da.SearchCursor)
    # =========================================================================

    def read_strike_dip_data(self, source, strike_field, dip_field,
                              overturned_field, label_field, strat_height_field,
                              dem_path, feedback):
        """Read strike, dip, and location data from a feature source."""
        has_z = QgsWkbTypes.hasZ(source.wkbType())

        # Open DEM for Z extraction
        dem_ds   = gdal.Open(dem_path)
        dem_gt   = dem_ds.GetGeoTransform() if dem_ds else None
        dem_band = dem_ds.GetRasterBand(1) if dem_ds else None
        dem_nd   = dem_band.GetNoDataValue() if dem_band else None

        data = {
            'x': [], 'y': [], 'z': [],
            'strike': [], 'dip': [],
            'overturned': [], 'labels': [],
            'strat_height': [],
        }

        for feat in source.getFeatures():
            geom = feat.geometry()
            pt   = geom.asPoint()
            data['x'].append(pt.x())
            data['y'].append(pt.y())

            # Z: from geometry or DEM
            z = None
            if has_z:
                v = geom.constGet()
                z = v.z() if hasattr(v, 'z') else None
                if z == 0.0:
                    z = None
            if z is None:
                z = self._sample_raster(dem_band, dem_gt, dem_nd, pt.x(), pt.y())
                if z is None:
                    feedback.pushWarning(
                        f'  Warning: Could not extract DEM value near ({pt.x():.1f}, {pt.y():.1f}), using 0.0'
                    )
                    z = 0.0
            data['z'].append(z)

            data['strike'].append(float(feat[strike_field]))
            data['dip'].append(float(feat[dip_field]))
            data['labels'].append(str(feat[label_field]))

            if overturned_field:
                val = feat[overturned_field]
                data['overturned'].append(1 if val == 1 else 0)
            else:
                data['overturned'].append(0)

            if strat_height_field:
                val = feat[strat_height_field]
                try:
                    data['strat_height'].append(float(val) if val is not None else None)
                except (TypeError, ValueError):
                    data['strat_height'].append(None)
            else:
                data['strat_height'].append(None)

        dem_ds = None

        for key in data:
            data[key] = np.array(data[key], dtype=object)

        return data

    # =========================================================================
    # calculate_fold_axis  (pure numpy — unchanged from arcpy version)
    # =========================================================================

    def calculate_fold_axis(self, strike_dip_data):
        """Calculate best-fit cylindrical fold axis from bedding measurements."""
        strikes = strike_dip_data['strike'].astype(float)
        dips    = strike_dip_data['dip'].astype(float)

        poles = []
        for strike, dip in zip(strikes, dips):
            pole_trend  = (strike - 90) % 360
            pole_plunge = 90 - dip

            pt_rad = np.radians(pole_trend)
            pp_rad = np.radians(pole_plunge)

            px = np.cos(pp_rad) * np.sin(pt_rad)
            py = np.cos(pp_rad) * np.cos(pt_rad)
            pz = np.sin(pp_rad)
            poles.append([px, py, pz])

        poles      = np.array(poles)
        covariance = np.dot(poles.T, poles)
        eigenvalues, eigenvectors = np.linalg.eig(covariance)

        min_idx  = np.argmin(eigenvalues)
        fold_axis = eigenvectors[:, min_idx]

        x, y, z = fold_axis
        plunge   = np.degrees(np.arcsin(float(z)))
        trend    = np.degrees(np.arctan2(float(x), float(y)))
        if trend < 0:
            trend += 360
        if plunge < 0:
            plunge = -plunge
            trend  = (trend + 180) % 360

        return trend, plunge

    # =========================================================================
    # calculate_projected_attitudes  (pure numpy — unchanged from arcpy version)
    # =========================================================================

    def calculate_projected_attitudes(self, strike_dip_data, fold_trend, fold_plunge):
        """Calculate the intersection of each bedding plane with the profile plane."""
        strikes    = strike_dip_data['strike'].astype(float)
        dips       = strike_dip_data['dip'].astype(float)
        overturned = strike_dip_data['overturned']

        trend_rad  = np.radians(fold_trend)
        plunge_rad = np.radians(fold_plunge)

        fold_axis = np.array([
            np.cos(plunge_rad) * np.sin(trend_rad),
            np.cos(plunge_rad) * np.cos(trend_rad),
            np.sin(plunge_rad),
        ])

        vertical  = np.array([0, 0, 1])
        v_par     = np.dot(vertical, fold_axis) * fold_axis
        profile_y = vertical - v_par
        profile_y = profile_y / np.linalg.norm(profile_y)

        profile_x = np.cross(profile_y, fold_axis)
        profile_x = profile_x / np.linalg.norm(profile_x)

        points_3d = np.column_stack([
            strike_dip_data['x'].astype(float),
            strike_dip_data['y'].astype(float),
            strike_dip_data['z'].astype(float),
        ])

        centroid        = np.mean(points_3d, axis=0)
        points_centered = points_3d - centroid

        profile_x_coords = np.dot(points_centered, profile_x)
        profile_y_coords = np.dot(points_centered, profile_y)

        projected_attitudes = []
        for i in range(len(strikes)):
            pole_trend  = strikes[i] + 90
            pole_plunge = 90 - dips[i]

            if overturned[i]:
                pole_trend  += 180
                pole_plunge  = -pole_plunge

            pole_trend = pole_trend % 360
            pt_rad = np.radians(pole_trend)
            pp_rad = np.radians(pole_plunge)

            pole = np.array([
                np.cos(pp_rad) * np.sin(pt_rad),
                np.cos(pp_rad) * np.cos(pt_rad),
                np.sin(pp_rad),
            ])

            intersection_3d = np.cross(pole, fold_axis)
            if np.linalg.norm(intersection_3d) > 1e-10:
                intersection_3d = intersection_3d / np.linalg.norm(intersection_3d)
            else:
                intersection_3d = profile_x.copy()

            ix = np.dot(intersection_3d, profile_x)
            iy = np.dot(intersection_3d, profile_y)

            norm = math.sqrt(ix ** 2 + iy ** 2)
            if norm > 1e-10:
                ix /= norm
                iy /= norm

            projected_attitudes.append({'x': ix, 'y': iy,
                                         'pole_3d': pole,
                                         'intersection_3d': intersection_3d})

        strike_dip_data['projected_attitudes'] = projected_attitudes
        strike_dip_data['profile_x']           = profile_x_coords
        strike_dip_data['profile_y']           = profile_y_coords

        # Stratigraphic height vectors
        strat_height_vectors = []
        for i in range(len(strikes)):
            att = projected_attitudes[i]
            is_overturned = overturned[i]

            perp_x = -att['y']
            perp_y =  att['x']
            norm   = math.sqrt(perp_x ** 2 + perp_y ** 2)
            if norm > 1e-10:
                perp_x /= norm
                perp_y /= norm

            if is_overturned:
                if perp_y > 0:
                    perp_x, perp_y = -perp_x, -perp_y
            else:
                if perp_y < 0:
                    perp_x, perp_y = -perp_x, -perp_y

            strat_height_vectors.append({'x': perp_x, 'y': perp_y})

        strike_dip_data['strat_height_vectors'] = strat_height_vectors
        return strike_dip_data

    # =========================================================================
    # analyze_wedges  (pure numpy — unchanged from arcpy version)
    # =========================================================================

    def analyze_wedges(self, strike_dip_data, sorted_indices):
        """Analyze wedge geometry between adjacent measurement pairs."""
        PARALLEL_TOLERANCE = 0.01
        wedge_data = []

        for i in range(len(sorted_indices) - 1):
            idx1 = sorted_indices[i]
            idx2 = sorted_indices[i + 1]

            wedge = {
                'point1_idx':        idx1,
                'point2_idx':        idx2,
                'type':              None,
                'intersection_point': None,
                'cylinder_axis_location': None,
            }

            att1 = strike_dip_data['projected_attitudes'][idx1]
            att2 = strike_dip_data['projected_attitudes'][idx2]

            dot_product = att1['x'] * att2['x'] + att1['y'] * att2['y']
            dot_product = np.clip(dot_product, -1.0, 1.0)
            angle_diff  = np.degrees(np.arccos(np.abs(dot_product)))

            if angle_diff < PARALLEL_TOLERANCE:
                sv1 = strike_dip_data['strat_height_vectors'][idx1]
                sv2 = strike_dip_data['strat_height_vectors'][idx2]
                if sv1['y'] * sv2['y'] >= 0:
                    wedge['type'] = 'parallel'
                else:
                    wedge['type'] = 'invalid_parallel'

            wedge_data.append(wedge)

        return wedge_data

    # =========================================================================
    # extract_dem_cell_centroids  (replaces arcpy.sa.ExtractByMask + RasterToNumPyArray)
    # =========================================================================

    def extract_dem_cell_centroids(self, dem_path, study_area_path, feedback):
        """Extract centroid coordinates and elevation values for DEM cells in the study area."""
        feedback.pushInfo('    Clipping DEM to study area...')
        ds, _ = self._open_clipped_raster(dem_path, study_area_path)
        if ds is None:
            feedback.pushWarning('Could not open clipped DEM')
            return None

        feedback.pushInfo('    Converting DEM to point array...')
        gt     = ds.GetGeoTransform()
        band   = ds.GetRasterBand(1)
        arr    = band.ReadAsArray()
        nodata = band.GetNoDataValue()
        if nodata is None:
            nodata = -9999

        nrows, ncols = arr.shape
        cell_x, cell_y, cell_z = [], [], []

        for row in range(nrows):
            for col in range(ncols):
                elev = float(arr[row, col])
                if elev == nodata:
                    continue
                # gt[3] = top-left Y, gt[5] < 0  →  centroid Y decreases with row
                cx = gt[0] + (col + 0.5) * gt[1]
                cy = gt[3] + (row + 0.5) * gt[5]
                cell_x.append(cx)
                cell_y.append(cy)
                cell_z.append(elev)

        ds = None
        return {
            'x': np.array(cell_x),
            'y': np.array(cell_y),
            'z': np.array(cell_z),
        }

    # =========================================================================
    # assign_dem_cells_to_wedges  (replaces arcpy SetProgressor; vectorized wedge case)
    # =========================================================================

    def assign_dem_cells_to_wedges(self, dem_profile_x, dem_profile_y,
                                    strike_dip_data, wedge_data, sorted_indices,
                                    feedback):
        """Assign each DEM cell to a wedge or rectangle."""
        feedback.pushInfo('  Assigning DEM cells to wedges/rectangles...')

        wedge_assignments = np.full(len(dem_profile_x), -1, dtype=int)
        valid_wedges      = [w for w in wedge_data if w['type'] in ('parallel', 'intersecting')]
        n_valid           = len(valid_wedges)

        for wedge_counter, wedge in enumerate(valid_wedges):
            feedback.setProgress(int(100 * wedge_counter / max(n_valid, 1)))

            left_idx  = wedge['point1_idx']
            right_idx = wedge['point2_idx']

            left_x  = strike_dip_data['profile_x'][left_idx]
            left_y  = strike_dip_data['profile_y'][left_idx]
            right_x = strike_dip_data['profile_x'][right_idx]
            right_y = strike_dip_data['profile_y'][right_idx]

            strat_vec_left = strike_dip_data['strat_height_vectors'][left_idx]
            wedge_idx      = wedge_data.index(wedge)

            if wedge['type'] == 'parallel':
                # RECTANGLE CASE — identical rotation logic to original
                dem_tx = dem_profile_x - left_x
                dem_ty = dem_profile_y - left_y
                rx     = right_x - left_x
                ry     = right_y - left_y

                angle = np.arctan2(strat_vec_left['x'], strat_vec_left['y'])
                cos_a, sin_a = np.cos(-angle), np.sin(-angle)

                dem_rot_x = cos_a * dem_tx - sin_a * dem_ty
                dem_rot_y = sin_a * dem_tx + cos_a * dem_ty
                right_rot_x = cos_a * rx - sin_a * ry
                right_rot_y = sin_a * rx + cos_a * ry

                in_x = ((dem_rot_x >= 0) & (dem_rot_x <= right_rot_x) if right_rot_x >= 0
                         else (dem_rot_x <= 0) & (dem_rot_x >= right_rot_x))
                in_y = ((dem_rot_y >= 0) & (dem_rot_y <= right_rot_y) if right_rot_y >= 0
                         else (dem_rot_y <= 0) & (dem_rot_y >= right_rot_y))
                in_region = in_x & in_y

            elif wedge['type'] == 'intersecting':
                # WEDGE CASE — vectorized version of original per-cell loop
                intersection = wedge['intersection_point']
                ix_pt, iy_pt = intersection[0], intersection[1]

                b1 = np.array([left_x - ix_pt, left_y - iy_pt])
                b2 = np.array([right_x - ix_pt, right_y - iy_pt])
                b1 /= np.linalg.norm(b1)
                b2 /= np.linalg.norm(b2)
                dot_between = float(np.dot(b1, b2))

                cvx = dem_profile_x - ix_pt
                cvy = dem_profile_y - iy_pt
                cv_len = np.sqrt(cvx ** 2 + cvy ** 2)
                nonzero = cv_len > 1e-10

                cvx_n = np.where(nonzero, cvx / (cv_len + 1e-30), 0.0)
                cvy_n = np.where(nonzero, cvy / (cv_len + 1e-30), 0.0)

                dot1 = cvx_n * b1[0] + cvy_n * b1[1]
                dot2 = cvx_n * b2[0] + cvy_n * b2[1]
                in_region = nonzero & (dot1 > dot_between) & (dot2 > dot_between)
            else:
                continue

            already   = (wedge_assignments >= 0) & in_region
            unassigned = (wedge_assignments == -1) & in_region
            wedge_assignments[already]   = -2
            wedge_assignments[unassigned] = wedge_idx

        feedback.setProgress(100)

        n_assigned = int(np.sum(wedge_assignments >= 0))
        n_multiple = int(np.sum(wedge_assignments == -2))
        feedback.pushInfo(f'  Assigned {n_assigned} DEM cells to wedges/rectangles')
        if n_multiple > 0:
            feedback.pushInfo(
                f'  Warning: {n_multiple} DEM cells belong to multiple wedges (marked as ambiguous)'
            )
        return wedge_assignments

    # =========================================================================
    # calculate_dem_stratigraphic_heights  (replaces arcpy SetProgressor)
    # =========================================================================

    def calculate_dem_stratigraphic_heights(self, dem_profile_x, dem_profile_y,
                                             dem_wedge_assignments, strike_dip_data,
                                             wedge_data, feedback):
        """Calculate stratigraphic heights for DEM cells based on their wedge assignments."""
        feedback.pushInfo('  Calculating stratigraphic heights for DEM cells...')

        dem_strat_heights = np.full(len(dem_profile_x), np.nan)

        if 'calculated_strat_height' not in strike_dip_data:
            feedback.pushWarning(
                '  Warning: No calculated stratigraphic heights available for measurements'
            )
            return dem_strat_heights

        point_heights = strike_dip_data['calculated_strat_height']
        valid_wedges  = [w for w in wedge_data if w['type'] in ('parallel', 'intersecting')]
        n_valid       = len(valid_wedges)

        for k, wedge in enumerate(valid_wedges):
            feedback.setProgress(int(100 * k / max(n_valid, 1)))
            wedge_idx = wedge_data.index(wedge)
            mask      = (dem_wedge_assignments == wedge_idx)
            if not np.any(mask):
                continue

            left_idx  = wedge['point1_idx']
            right_idx = wedge['point2_idx']

            left_height = point_heights[left_idx]
            if left_height is None:
                feedback.pushWarning(
                    f'  Warning: No height available for left point in wedge {wedge_idx}'
                )
                continue

            left_x = strike_dip_data['profile_x'][left_idx]
            left_y = strike_dip_data['profile_y'][left_idx]

            if wedge['type'] == 'parallel':
                strat_vec = strike_dip_data['strat_height_vectors'][left_idx]
                dem_x     = dem_profile_x[mask]
                dem_y     = dem_profile_y[mask]

                dem_tx = dem_x - left_x
                dem_ty = dem_y - left_y

                angle = np.arctan2(strat_vec['x'], strat_vec['y'])
                cos_a, sin_a = np.cos(-angle), np.sin(-angle)

                dem_rot_y = sin_a * dem_tx + cos_a * dem_ty
                dem_strat_heights[mask] = dem_rot_y + left_height

            elif wedge['type'] == 'intersecting':
                intersection      = wedge['intersection_point']
                intersection_type = wedge.get('intersection_type', 'unknown')

                dem_x = dem_profile_x[mask]
                dem_y = dem_profile_y[mask]

                dx_dem = dem_x - intersection[0]
                dy_dem = dem_y - intersection[1]
                dist_from_int = np.sqrt(dx_dem ** 2 + dy_dem ** 2)

                dx_left = left_x - intersection[0]
                dy_left = left_y - intersection[1]
                dist_left = math.sqrt(dx_left ** 2 + dy_left ** 2)

                if intersection_type == 'up-up':
                    constant = left_height + dist_left
                    dem_strat_heights[mask] = -dist_from_int + constant
                elif intersection_type == 'down-down':
                    constant = left_height - dist_left
                    dem_strat_heights[mask] = dist_from_int + constant
                else:
                    feedback.pushWarning(
                        f'  Warning: Unknown intersection type for wedge {wedge_idx}'
                    )

        feedback.setProgress(100)
        n_calc = int(np.sum(~np.isnan(dem_strat_heights)))
        feedback.pushInfo(f'  Calculated stratigraphic heights for {n_calc} DEM cells')
        return dem_strat_heights

    # =========================================================================
    # create_wedge_assignment_raster  (replaces arcpy.NumPyArrayToRaster)
    # =========================================================================

    def create_wedge_assignment_raster(self, dem_cell_data, wedge_assignments,
                                        dem_raster, study_area, output_path, feedback):
        """Create a raster showing which wedge each DEM cell belongs to."""
        ds, _ = self._open_clipped_raster(dem_raster, study_area)
        gt    = ds.GetGeoTransform()
        proj  = ds.GetProjection()
        band  = ds.GetRasterBand(1)
        nrows = ds.RasterYSize
        ncols = ds.RasterXSize
        ds    = None

        wedge_array = np.full((nrows, ncols), -9999, dtype=np.int16)

        for i in range(len(dem_cell_data['x'])):
            x   = dem_cell_data['x'][i]
            y   = dem_cell_data['y'][i]
            wid = int(wedge_assignments[i])

            col = int((x - gt[0]) / gt[1])
            row = int((y - gt[3]) / gt[5])   # gt[5] < 0 → correct

            if 0 <= row < nrows and 0 <= col < ncols:
                wedge_array[row, col] = wid

        drv    = gdal.GetDriverByName('GTiff')
        out_ds = drv.Create(output_path, ncols, nrows, 1, gdal.GDT_Int16)
        out_ds.SetGeoTransform(gt)
        out_ds.SetProjection(proj)
        out_band = out_ds.GetRasterBand(1)
        out_band.WriteArray(wedge_array)
        out_band.SetNoDataValue(-9999)
        out_ds.FlushCache()
        out_ds = None

    # =========================================================================
    # calculate_stratigraphic_heights_in_profile  (replaces arcpy AddMessage/Warning)
    # =========================================================================

    def calculate_stratigraphic_heights_in_profile(self, strike_dip_data, wedge_data,
                                                    profile_x_coords, profile_y_coords,
                                                    sorted_indices, feedback):
        """Calculate stratigraphic heights for measurements using the Busk method."""
        strat_heights = strike_dip_data['strat_height']
        known_heights = [h for h in strat_heights if h is not None]

        if len(known_heights) == 0:
            feedback.pushWarning('No stratigraphic height specified - cannot calculate heights')
            return
        if len(known_heights) > 1:
            feedback.pushWarning(
                f'Multiple measurements have stratigraphic heights specified ({len(known_heights)} found)'
            )
            feedback.pushWarning(
                'Only one measurement should have a known stratigraphic height - cannot calculate heights'
            )
            return

        known_idx, known_height = None, None
        for i, h in enumerate(strat_heights):
            if h is not None:
                known_idx    = i
                known_height = float(h)
                break

        feedback.pushInfo(
            f"  Using {strike_dip_data['labels'][known_idx]} as reference "
            f'with height = {known_height:.2f} m'
        )

        calculated_heights            = [None] * len(strike_dip_data['x'])
        leftmost_idx                  = sorted_indices[0]
        calculated_heights[leftmost_idx] = 0.0
        feedback.pushInfo(
            f"  Starting propagation from leftmost point "
            f"{strike_dip_data['labels'][leftmost_idx]} with initial height = 0.0"
        )

        for i in range(len(sorted_indices) - 1):
            left_idx  = sorted_indices[i]
            right_idx = sorted_indices[i + 1]

            wedge = None
            for w in wedge_data:
                if ((w['point1_idx'] == left_idx and w['point2_idx'] == right_idx) or
                        (w['point1_idx'] == right_idx and w['point2_idx'] == left_idx)):
                    wedge = w
                    break

            if wedge is None or wedge['type'] not in ('parallel', 'intersecting'):
                feedback.pushWarning(
                    f"  No valid wedge found between "
                    f"{strike_dip_data['labels'][left_idx]} and "
                    f"{strike_dip_data['labels'][right_idx]}"
                )
                continue

            left_height = calculated_heights[left_idx]
            left_x      = float(profile_x_coords[left_idx])
            left_y      = float(profile_y_coords[left_idx])
            right_x     = float(profile_x_coords[right_idx])
            right_y     = float(profile_y_coords[right_idx])

            if wedge['type'] == 'intersecting':
                ix, iy            = wedge['intersection_point']
                intersection_type = wedge.get('intersection_type', None)

                dist_left  = math.sqrt((left_x - ix)  ** 2 + (left_y - iy)  ** 2)
                dist_right = math.sqrt((right_x - ix) ** 2 + (right_y - iy) ** 2)

                feedback.pushInfo(f'    Intersection at ({ix:.2f}, {iy:.2f})')
                feedback.pushInfo(f'    Left point ({left_x:.2f}, {left_y:.2f}), dist={dist_left:.2f}')
                feedback.pushInfo(f'    Right point ({right_x:.2f}, {right_y:.2f}), dist={dist_right:.2f}')

                if intersection_type == 'up-up':
                    constant     = left_height + dist_left
                    right_height = -dist_right + constant
                    feedback.pushInfo(
                        f"    UP-UP wedge: {strike_dip_data['labels'][left_idx]} ({left_height:.2f})"
                        f" -> {strike_dip_data['labels'][right_idx]} ({right_height:.2f})"
                    )
                elif intersection_type == 'down-down':
                    constant     = left_height - dist_left
                    right_height = dist_right + constant
                    feedback.pushInfo(
                        f"    DOWN-DOWN wedge: {strike_dip_data['labels'][left_idx]} ({left_height:.2f})"
                        f" -> {strike_dip_data['labels'][right_idx]} ({right_height:.2f})"
                    )
                else:
                    feedback.pushWarning(
                        f"    Unknown intersection type for wedge between "
                        f"{strike_dip_data['labels'][left_idx]} and "
                        f"{strike_dip_data['labels'][right_idx]}"
                    )
                    continue

                calculated_heights[right_idx] = right_height

            elif wedge['type'] == 'parallel':
                sv  = strike_dip_data['strat_height_vectors'][left_idx]
                dx  = right_x - left_x
                dy  = right_y - left_y
                proj = dx * sv['x'] + dy * sv['y']
                right_height = left_height + proj

                direction = 'up' if proj > 0 else 'down'
                feedback.pushInfo(
                    f"    PARALLEL ({direction}): "
                    f"{strike_dip_data['labels'][left_idx]} ({left_height:.2f})"
                    f" -> {strike_dip_data['labels'][right_idx]} ({right_height:.2f}), "
                    f"offset={proj:.2f}"
                )
                calculated_heights[right_idx] = right_height

        if calculated_heights[known_idx] is not None:
            correction = known_height - calculated_heights[known_idx]
            feedback.pushInfo(f'\n  Applying correction constant: {correction:.2f} m')
            for i in range(len(calculated_heights)):
                if calculated_heights[i] is not None:
                    calculated_heights[i] += correction
            feedback.pushInfo('\n  Final stratigraphic heights:')
            for i in range(len(calculated_heights)):
                if calculated_heights[i] is not None:
                    feedback.pushInfo(
                        f"    {strike_dip_data['labels'][i]}: {calculated_heights[i]:.2f} m"
                    )
        else:
            feedback.pushWarning(
                f"  Could not propagate to reference point "
                f"{strike_dip_data['labels'][known_idx]}"
            )

        strike_dip_data['calculated_strat_height'] = calculated_heights

    # =========================================================================
    # create_output_strike_dip_fc  (replaces arcpy feature class I/O)
    # =========================================================================

    def create_output_strike_dip_fc(self, source, strike_dip_data, output_sink,
                                     label_field, feedback):
        """Write output feature class with calculated stratigraphic heights appended."""
        if output_sink is None:
            return

        feedback.pushInfo(f'  Creating output feature class...')

        calculated_heights = strike_dip_data.get('calculated_strat_height', [])
        labels             = strike_dip_data['labels']

        label_to_height = {}
        for i, lbl in enumerate(labels):
            if i < len(calculated_heights) and calculated_heights[i] is not None:
                label_to_height[str(lbl)] = calculated_heights[i]

        feedback.pushInfo(f'  Created mapping for {len(label_to_height)} features with heights')

        update_count = 0
        for feat in source.getFeatures():
            new_feat = QgsFeature(feat)
            attrs    = list(feat.attributes())
            lbl_val  = feat[label_field]
            ht       = label_to_height.get(str(lbl_val), None)
            attrs.append(ht)
            new_feat.setAttributes(attrs)
            output_sink.addFeature(new_feat, QgsFeatureSink.FastInsert)
            if ht is not None:
                update_count += 1

        feedback.pushInfo(f'  Updated {update_count} features with calculated heights')

    # =========================================================================
    # create_strat_height_raster  (replaces arcpy.NumPyArrayToRaster)
    # =========================================================================

    def create_strat_height_raster(self, dem_cell_data, dem_strat_heights,
                                    dem_raster, study_area, output_path, feedback):
        """Create a raster showing stratigraphic height for each DEM cell."""
        ds, _ = self._open_clipped_raster(dem_raster, study_area)
        gt    = ds.GetGeoTransform()
        proj  = ds.GetProjection()
        nrows = ds.RasterYSize
        ncols = ds.RasterXSize
        ds    = None

        strat_array = np.full((nrows, ncols), -9999.0, dtype=np.float32)

        for i in range(len(dem_cell_data['x'])):
            if np.isnan(dem_strat_heights[i]):
                continue
            x = dem_cell_data['x'][i]
            y = dem_cell_data['y'][i]

            col = int((x - gt[0]) / gt[1])
            row = int((y - gt[3]) / gt[5])

            if 0 <= row < nrows and 0 <= col < ncols:
                strat_array[row, col] = float(dem_strat_heights[i])

        drv    = gdal.GetDriverByName('GTiff')
        out_ds = drv.Create(output_path, ncols, nrows, 1, gdal.GDT_Float32)
        out_ds.SetGeoTransform(gt)
        out_ds.SetProjection(proj)
        out_band = out_ds.GetRasterBand(1)
        out_band.WriteArray(strat_array)
        out_band.SetNoDataValue(-9999.0)
        out_ds.FlushCache()
        out_ds = None

    # =========================================================================
    # create_downplunge_view  (replaces arcpy orchestration + matplotlib output)
    # =========================================================================

    def create_downplunge_view(self, strike_dip_data, wedge_data, dem_cell_data,
                                fold_trend, fold_plunge,
                                output_pdf, dem_raster, study_area,
                                wedge_raster_output,
                                strike_dip_source, output_sink, label_field,
                                plot_dem_cells, strat_height_raster_output,
                                feedback):
        """Create downplunge projection, run analysis, and save as PDF."""

        # Re-derive profile plane basis vectors (same as calculate_projected_attitudes)
        trend_rad  = np.radians(fold_trend)
        plunge_rad = np.radians(fold_plunge)
        fold_axis  = np.array([
            np.cos(plunge_rad) * np.sin(trend_rad),
            np.cos(plunge_rad) * np.cos(trend_rad),
            -np.sin(plunge_rad),   # sign matches create_downplunge_view in original
        ])

        vertical  = np.array([0, 0, 1])
        v_par     = np.dot(vertical, fold_axis) * fold_axis
        profile_y = vertical - v_par
        profile_y = profile_y / np.linalg.norm(profile_y)
        profile_x = np.cross(profile_y, fold_axis)
        profile_x = profile_x / np.linalg.norm(profile_x)

        # Re-project measurement points (centred)
        points_3d = np.column_stack([
            strike_dip_data['x'].astype(float),
            strike_dip_data['y'].astype(float),
            strike_dip_data['z'].astype(float),
        ])
        centroid        = np.mean(points_3d, axis=0)
        points_centered = points_3d - centroid

        profile_x_coords = np.dot(points_centered, profile_x)
        profile_y_coords = np.dot(points_centered, profile_y)

        # Update stored profile coordinates
        strike_dip_data['profile_x'] = profile_x_coords
        strike_dip_data['profile_y'] = profile_y_coords

        sorted_indices = np.argsort(profile_x_coords)

        # Project DEM cells onto profile plane
        dem_profile_x = None
        dem_profile_y = None
        if dem_cell_data is not None:
            feedback.pushInfo('  Projecting DEM cells onto profile plane...')
            dem_pts = np.column_stack([
                dem_cell_data['x'],
                dem_cell_data['y'],
                dem_cell_data['z'],
            ])
            dem_pts_centered = dem_pts - centroid
            dem_profile_x    = np.dot(dem_pts_centered, profile_x)
            dem_profile_y    = np.dot(dem_pts_centered, profile_y)

        # Calculate wedge intersections and update wedge_data in place
        x_range = float(profile_x_coords.max() - profile_x_coords.min())

        for wedge in wedge_data:
            if wedge['type'] == 'parallel':
                continue

            if wedge['type'] == 'invalid_parallel':
                idx1, idx2 = wedge['point1_idx'], wedge['point2_idx']
                feedback.pushWarning(
                    f"No intersection between stratigraphic up or down for measurements "
                    f"{strike_dip_data['labels'][idx1]} and {strike_dip_data['labels'][idx2]}"
                )
                continue

            idx1, idx2 = wedge['point1_idx'], wedge['point2_idx']
            label1     = strike_dip_data['labels'][idx1]
            label2     = strike_dip_data['labels'][idx2]

            sp1 = int(np.where(sorted_indices == idx1)[0][0])
            sp2 = int(np.where(sorted_indices == idx2)[0][0])
            if abs(sp1 - sp2) != 1:
                continue

            p1_x = float(profile_x_coords[idx1])
            p1_y = float(profile_y_coords[idx1])
            p2_x = float(profile_x_coords[idx2])
            p2_y = float(profile_y_coords[idx2])

            sv1 = strike_dip_data['strat_height_vectors'][idx1]
            sv2 = strike_dip_data['strat_height_vectors'][idx2]

            feedback.pushInfo(f'\nChecking intersection for {label1} and {label2}:')
            feedback.pushInfo(f"  {label1} strat height vector: x={sv1['x']:+.4f}, y={sv1['y']:+.4f}")
            feedback.pushInfo(f"  {label2} strat height vector: x={sv2['x']:+.4f}, y={sv2['y']:+.4f}")

            intersection, t1, t2 = self.find_line_intersection_with_params(
                p1_x, p1_y, sv1['x'], sv1['y'],
                p2_x, p2_y, sv2['x'], sv2['y'],
            )

            if intersection is None:
                feedback.pushInfo('  -> No intersection (parallel vectors)')
                feedback.pushWarning(
                    f'No intersection between stratigraphic height vectors for '
                    f'measurements {label1} and {label2}'
                )
                continue

            feedback.pushInfo(f'  -> Intersection at ({intersection[0]:.2f}, {intersection[1]:.2f})')
            feedback.pushInfo(f'  -> Parameters: t1={t1:.2f}, t2={t2:.2f}')

            if t1 > 0 and t2 > 0:
                wedge['type']               = 'intersecting'
                wedge['intersection_point'] = intersection
                wedge['cylinder_axis_location'] = intersection
                wedge['intersection_type']  = 'up-up'
                feedback.pushInfo('  -> UP-UP wedge (both t positive)')
            elif t1 < 0 and t2 < 0:
                wedge['type']               = 'intersecting'
                wedge['intersection_point'] = intersection
                wedge['cylinder_axis_location'] = intersection
                wedge['intersection_type']  = 'down-down'
                feedback.pushInfo('  -> DOWN-DOWN wedge (both t negative)')
            else:
                feedback.pushInfo(
                    f'  -> INVALID: t values have opposite signs (t1={t1:.2f}, t2={t2:.2f})'
                )
                feedback.pushWarning(
                    f'No valid wedge geometry for measurements {label1} and {label2}'
                )

        # Validate coverage
        feedback.pushInfo('\nValidating wedge coverage...')
        all_valid = all(w['type'] in ('parallel', 'intersecting') for w in wedge_data)

        if not all_valid:
            feedback.pushWarning(
                'Stratigraphic heights cannot be calculated in the study area because '
                'some areas are not covered by valid wedges or rectangles.'
            )
            feedback.pushWarning(
                'All adjacent measurement pairs must have either valid intersections '
                'or parallel geometry with consistent stratigraphic directions.'
            )
        else:
            feedback.pushInfo('  All adjacent pairs have valid geometry ✓')

            # Calculate stratigraphic heights
            if any(h is not None for h in strike_dip_data['strat_height']):
                feedback.pushInfo('\nCalculating stratigraphic heights...')
                self.calculate_stratigraphic_heights_in_profile(
                    strike_dip_data, wedge_data,
                    profile_x_coords, profile_y_coords, sorted_indices, feedback,
                )
            else:
                feedback.pushInfo(
                    '\nNo stratigraphic height field provided - '
                    'assuming leftmost point has stratigraphic height = 0'
                )
                n_pts         = len(strike_dip_data['x'])
                strat_height  = [None] * n_pts
                strat_height[sorted_indices[0]] = 0.0
                strike_dip_data['strat_height'] = np.array(strat_height, dtype=object)
                feedback.pushInfo('\nCalculating stratigraphic heights...')
                self.calculate_stratigraphic_heights_in_profile(
                    strike_dip_data, wedge_data,
                    profile_x_coords, profile_y_coords, sorted_indices, feedback,
                )

            # Output feature class
            if output_sink is not None:
                feedback.pushInfo('\nCreating output strike/dip feature class with calculated heights...')
                self.create_output_strike_dip_fc(
                    strike_dip_source, strike_dip_data, output_sink, label_field, feedback,
                )

            # Assign DEM cells to wedges and build rasters
            dem_wedge_assignments = None
            dem_strat_heights     = None
            if dem_profile_x is not None:
                dem_wedge_assignments = self.assign_dem_cells_to_wedges(
                    dem_profile_x, dem_profile_y,
                    strike_dip_data, wedge_data, sorted_indices, feedback,
                )
                dem_strat_heights = self.calculate_dem_stratigraphic_heights(
                    dem_profile_x, dem_profile_y,
                    dem_wedge_assignments, strike_dip_data, wedge_data, feedback,
                )

                if wedge_raster_output:
                    feedback.pushInfo('\nCreating wedge assignment raster...')
                    self.create_wedge_assignment_raster(
                        dem_cell_data, dem_wedge_assignments,
                        dem_raster, study_area, wedge_raster_output, feedback,
                    )
                    feedback.pushInfo(f'  Wedge assignment raster saved to: {wedge_raster_output}')

                if strat_height_raster_output:
                    feedback.pushInfo('\nCreating stratigraphic height raster...')
                    self.create_strat_height_raster(
                        dem_cell_data, dem_strat_heights,
                        dem_raster, study_area, strat_height_raster_output, feedback,
                    )
                    feedback.pushInfo(
                        f'  Stratigraphic height raster saved to: {strat_height_raster_output}'
                    )

        # ── Build the plot ───────────────────────────────────────────────────
        x_range = float(profile_x_coords.max() - profile_x_coords.min())
        y_range = float(profile_y_coords.max() - profile_y_coords.min())

        # Filter intersections that are too far away (> 5× plot width)
        max_distance = 5 * x_range
        center_x     = (profile_x_coords.min() + profile_x_coords.max()) / 2
        center_y     = (profile_y_coords.min() + profile_y_coords.max()) / 2
        for wedge in wedge_data:
            if wedge.get('intersection_point') is not None:
                ix, iy = wedge['intersection_point']
                dist   = math.sqrt((ix - center_x) ** 2 + (iy - center_y) ** 2)
                wedge['plot_intersection'] = (dist <= max_distance)

        fig, ax = plt.subplots(figsize=(10, 8))

        # Plot DEM cells coloured by wedge assignment
        if (dem_profile_x is not None and dem_profile_y is not None and plot_dem_cells
                and 'dem_wedge_assignments' in dir()):
            n_dem = len(dem_profile_x)
            feedback.pushInfo(f'  Plotting DEM cells...')
            max_plot = 10000
            if n_dem > max_plot:
                feedback.pushInfo(f'  Subsampling {max_plot} of {n_dem} DEM cells for plotting...')
                idx_plot  = np.random.choice(n_dem, max_plot, replace=False)
            else:
                feedback.pushInfo(f'  Plotting all {n_dem} DEM cells...')
                idx_plot = np.arange(n_dem)

            dem_px   = dem_profile_x[idx_plot]
            dem_py   = dem_profile_y[idx_plot]
            dem_wa   = dem_wedge_assignments[idx_plot]
            n_wedges = len(wedge_data)
            cmap     = plt.cm.get_cmap('tab20', n_wedges)

            for wi in range(n_wedges):
                m = (dem_wa == wi)
                if m.any():
                    ax.scatter(dem_px[m], dem_py[m], c=[cmap(wi)],
                               s=0.5, alpha=0.5, zorder=1, label=f'Wedge {wi}')

            m_un = (dem_wa == -1)
            if m_un.any():
                ax.scatter(dem_px[m_un], dem_py[m_un],
                           c='lightgray', s=0.5, alpha=0.3, zorder=1, label='Unassigned')
            m_am = (dem_wa == -2)
            if m_am.any():
                ax.scatter(dem_px[m_am], dem_py[m_am],
                           c='black', s=1.0, alpha=0.7, zorder=1,
                           label='Ambiguous (multiple wedges)')

        line_length = 0.03 * x_range

        # Measurement points
        ax.scatter(profile_x_coords, profile_y_coords,
                   c='red', s=25, edgecolors='black', linewidths=0.5,
                   zorder=3, label='Measurements')

        # Bedding traces and strat height vectors
        for i in range(len(profile_x_coords)):
            att = strike_dip_data['projected_attitudes'][i]
            sv  = strike_dip_data['strat_height_vectors'][i]
            dx  = att['x'] * line_length / 2
            dy  = att['y'] * line_length / 2

            ax.plot([profile_x_coords[i] - dx, profile_x_coords[i] + dx],
                    [profile_y_coords[i] - dy, profile_y_coords[i] + dy],
                    'b-', linewidth=1.5, zorder=2,
                    label='Bedding Trace' if i == 0 else '')

            ax.plot([profile_x_coords[i],
                     profile_x_coords[i] + sv['x'] * line_length / 2],
                    [profile_y_coords[i],
                     profile_y_coords[i] + sv['y'] * line_length / 2],
                    color='gold', linewidth=1.5, zorder=2,
                    label='Strat Height Vector' if i == 0 else '')

        # Labels
        for i in range(len(profile_x_coords)):
            ax.annotate(
                str(strike_dip_data['labels'][i]),
                (profile_x_coords[i], profile_y_coords[i]),
                xytext=(3, 3), textcoords='offset points', fontsize=8,
            )

        # Wedge geometry
        for wedge in wedge_data:
            idx1, idx2 = wedge['point1_idx'], wedge['point2_idx']
            sp1 = int(np.where(sorted_indices == idx1)[0][0])
            sp2 = int(np.where(sorted_indices == idx2)[0][0])
            if abs(sp1 - sp2) != 1:
                continue

            p1x = float(profile_x_coords[idx1])
            p1y = float(profile_y_coords[idx1])
            p2x = float(profile_x_coords[idx2])
            p2y = float(profile_y_coords[idx2])

            if wedge['type'] == 'parallel':
                att1 = strike_dip_data['projected_attitudes'][idx1]
                dx   = att1['x'] * line_length / 2
                dy   = att1['y'] * line_length / 2
                cx   = [p1x - dx, p1x + dx, p2x + dx, p2x - dx, p1x - dx]
                cy   = [p1y - dy, p1y + dy, p2y + dy, p2y - dy, p1y - dy]
                ax.plot(cx, cy, 'g-', linewidth=1, alpha=0.7, zorder=1)

            elif wedge['type'] == 'intersecting' and wedge.get('plot_intersection', False):
                ix, iy = wedge['intersection_point']
                ax.plot([p1x, ix], [p1y, iy], 'k--', linewidth=0.8, alpha=0.6, zorder=1)
                ax.plot([p2x, ix], [p2y, iy], 'k--', linewidth=0.8, alpha=0.6, zorder=1)
                ax.plot(ix, iy, 'ko', markersize=4, zorder=2)

        # Axes formatting
        padding = 0.1
        ax.set_xlim(profile_x_coords.min() - padding * x_range,
                    profile_x_coords.max() + padding * x_range)
        ax.set_ylim(profile_y_coords.min() - padding * y_range,
                    profile_y_coords.max() + padding * y_range)
        ax.set_xlabel('Profile Distance (m)', fontsize=12)
        ax.set_ylabel('Vertical Distance (m)', fontsize=12)
        ax.set_title('Downplunge Projection of Strike/Dip Measurements',
                     fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.set_aspect('equal', adjustable='datalim')
        ax.legend(loc='best')
        ax.invert_xaxis()

        self.add_scale_bar(ax, x_range)

        ax.text(0.02, 0.98,
                f'Fold Axis: {fold_trend:.1f}° / {fold_plunge:.1f}°',
                transform=ax.transAxes, fontsize=11, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

        with PdfPages(output_pdf) as pdf:
            pdf.savefig(fig, bbox_inches='tight')
            plt.close()

    # =========================================================================
    # Geometric utilities  (pure Python — identical to arcpy version)
    # =========================================================================

    def find_line_intersection(self, p1x, p1y, d1x, d1y, p2x, p2y, d2x, d2y):
        """Find intersection of two 2D lines; returns (x, y) or None if parallel."""
        det = d1x * (-d2y) - (-d2x) * d1y
        if abs(det) < 1e-10:
            return None
        dx = p2x - p1x
        dy = p2y - p1y
        t1 = (dx * (-d2y) - (-d2x) * dy) / det
        return (p1x + t1 * d1x, p1y + t1 * d1y)

    def find_line_intersection_with_params(self, p1x, p1y, d1x, d1y,
                                            p2x, p2y, d2x, d2y):
        """2D line intersection; returns ((x, y), t1, t2) or (None, None, None)."""
        det = d1x * (-d2y) - (-d2x) * d1y
        if abs(det) < 1e-10:
            return None, None, None
        dx = p2x - p1x
        dy = p2y - p1y
        t1 = (dx * (-d2y) - (-d2x) * dy) / det
        t2 = (d1x * dy - dx * d1y) / det
        return (p1x + t1 * d1x, p1y + t1 * d1y), t1, t2

    # =========================================================================
    # add_scale_bar  (pure matplotlib — identical to arcpy version)
    # =========================================================================

    def add_scale_bar(self, ax, x_range):
        """Add a scale bar to the lower-right corner of the plot."""
        nice_numbers  = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000]
        target_length = x_range * 0.25
        scale_length  = min(nice_numbers, key=lambda v: abs(v - target_length))

        xlim    = ax.get_xlim()
        ylim    = ax.get_ylim()
        y_range = ylim[1] - ylim[0]

        bar_x_start = xlim[1] + 0.05 * x_range
        bar_x_end   = bar_x_start + scale_length
        bar_y       = ylim[0] + 0.08 * y_range

        ax.plot([bar_x_start, bar_x_end], [bar_y, bar_y], 'k-', linewidth=2, zorder=10)

        tick_h = 0.01 * y_range
        for bx in (bar_x_start, bar_x_end):
            ax.plot([bx, bx], [bar_y - tick_h, bar_y + tick_h], 'k-', linewidth=2, zorder=10)

        ax.text(
            (bar_x_start + bar_x_end) / 2,
            bar_y + 3 * tick_h,
            f'{scale_length} m',
            ha='center', va='bottom', fontsize=10, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8),
        )

    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

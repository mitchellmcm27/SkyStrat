"""
SkyStrat Toolbox for Stratigraphic Height Calculation
Calculates stratigraphic heights from strike/dip measurements using the Busk cross-section method
"""

import arcpy
import numpy as np
import os
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


class Toolbox(object):
    def __init__(self):
        """Define the toolbox"""
        self.label = "SkyStrat Toolbox"
        self.alias = "skystrat"
        self.tools = [SkyStratDownplungeView]


class SkyStratDownplungeView(object):
    def __init__(self):
        """Define the tool"""
        self.label = "SkyStrat Downplunge View"
        self.description = "Creates a downplunge projection of strike/dip measurements"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        
        # Parameter 0: Input strike/dip point feature class
        param0 = arcpy.Parameter(
            displayName="Strike and Dip Point Feature Class",
            name="strike_dip_fc",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
        param0.filter.list = ["Point"]
        
        # Parameter 1: Strike field
        param1 = arcpy.Parameter(
            displayName="Strike Field",
            name="strike_field",
            datatype="Field",
            parameterType="Required",
            direction="Input")
        param1.parameterDependencies = [param0.name]
        
        # Parameter 2: Dip field
        param2 = arcpy.Parameter(
            displayName="Dip Field",
            name="dip_field",
            datatype="Field",
            parameterType="Required",
            direction="Input")
        param2.parameterDependencies = [param0.name]
        
        # Parameter 3: Overturned flag field (optional)
        param3 = arcpy.Parameter(
            displayName="Overturned Flag Field (Optional)",
            name="overturned_field",
            datatype="Field",
            parameterType="Optional",
            direction="Input")
        param3.parameterDependencies = [param0.name]
        
        # Parameter 4: Label field
        param4 = arcpy.Parameter(
            displayName="Label Field (for point identification)",
            name="label_field",
            datatype="Field",
            parameterType="Required",
            direction="Input")
        param4.parameterDependencies = [param0.name]
        param4.value = "OBJECTID"
        
        # Parameter 5: Stratigraphic height field (optional)
        param5 = arcpy.Parameter(
            displayName="Stratigraphic Height Field (Optional)",
            name="strat_height_field",
            datatype="Field",
            parameterType="Optional",
            direction="Input")
        param5.parameterDependencies = [param0.name]
        
        # Parameter 6: DEM raster
        param6 = arcpy.Parameter(
            displayName="Digital Elevation Model (DEM)",
            name="dem_raster",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input")
        
        # Parameter 7: Study area polygon
        param7 = arcpy.Parameter(
            displayName="Study Area Polygon",
            name="study_area",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
        param7.filter.list = ["Polygon"]
        
        # Parameter 8: Fold axis trend (optional)
        param8 = arcpy.Parameter(
            displayName="Fold Axis Trend (azimuth, optional)",
            name="fold_trend",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")
        
        # Parameter 9: Fold axis plunge (optional)
        param9 = arcpy.Parameter(
            displayName="Fold Axis Plunge (degrees, optional)",
            name="fold_plunge",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")
        
        # Parameter 10: Output PDF path
        param10 = arcpy.Parameter(
            displayName="Output PDF Path",
            name="output_pdf",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        param10.filter.list = ["pdf"]
        
        # Parameter 11: Plot DEM cell centroids
        param11 = arcpy.Parameter(
            displayName="Plot DEM Cell Centroids on Profile",
            name="plot_dem_cells",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")
        param11.value = False
        
        # Parameter 12: Wedge assignment raster output (optional)
        param12 = arcpy.Parameter(
            displayName="Output Wedge Assignment Raster (Optional)",
            name="wedge_raster_output",
            datatype="DERasterDataset",
            parameterType="Optional",
            direction="Output")
        
        # Parameter 13: Output strike/dip feature class with calculated heights (optional)
        param13 = arcpy.Parameter(
            displayName="Output Strike/Dip Feature Class with Heights (Optional)",
            name="output_strike_dip_fc",
            datatype="DEFeatureClass",
            parameterType="Optional",
            direction="Output")
        
        # Parameter 14: Stratigraphic height raster output (optional)
        param14 = arcpy.Parameter(
            displayName="Output Stratigraphic Height Raster (Optional)",
            name="strat_height_raster_output",
            datatype="DERasterDataset",
            parameterType="Optional",
            direction="Output")
        
        params = [param0, param1, param2, param3, param4, param5, param6, param7, param8, param9, param10, param11, param12, param13, param14]
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal validation"""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation"""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        
        # Get parameters
        strike_dip_fc = parameters[0].valueAsText
        strike_field = parameters[1].valueAsText
        dip_field = parameters[2].valueAsText
        overturned_field = parameters[3].valueAsText
        label_field = parameters[4].valueAsText
        strat_height_field = parameters[5].valueAsText
        dem_raster = parameters[6].valueAsText
        study_area = parameters[7].valueAsText
        fold_trend = parameters[8].value
        fold_plunge = parameters[9].value
        output_pdf = parameters[10].valueAsText
        plot_dem_cells = parameters[11].value
        wedge_raster_output = parameters[12].valueAsText
        output_strike_dip_fc = parameters[13].valueAsText
        strat_height_raster_output = parameters[14].valueAsText
        
        arcpy.AddMessage("=" * 60)
        arcpy.AddMessage("SkyStrat - Downplunge View Generation")
        arcpy.AddMessage("=" * 60)
        
        # Check coordinate systems
        self.check_utm_projection(strike_dip_fc, "Strike/Dip Feature Class")
        self.check_utm_projection(dem_raster, "DEM Raster")
        self.check_utm_projection(study_area, "Study Area Polygon")
        
        # Read strike and dip data
        arcpy.AddMessage("\nReading strike and dip measurements...")
        strike_dip_data = self.read_strike_dip_data(
            strike_dip_fc, strike_field, dip_field, overturned_field, label_field, strat_height_field, dem_raster)
        arcpy.AddMessage(f"  Loaded {len(strike_dip_data['x'])} measurements")
        
        # Calculate or use provided fold axis
        if fold_trend is None or fold_plunge is None:
            arcpy.AddMessage("\nCalculating best-fit fold axis...")
            fold_trend, fold_plunge = self.calculate_fold_axis(strike_dip_data)
            arcpy.AddMessage(f"  Calculated fold axis: {fold_trend:.1f}° / {fold_plunge:.1f}°")
        else:
            arcpy.AddMessage(f"\nUsing provided fold axis: {fold_trend:.1f}° / {fold_plunge:.1f}°")
        
        # Calculate projected attitudes on profile plane
        arcpy.AddMessage("\nCalculating projected attitudes on profile plane...")
        strike_dip_data = self.calculate_projected_attitudes(strike_dip_data, fold_trend, fold_plunge)
        arcpy.AddMessage(f"  Projected {len(strike_dip_data['x'])} attitudes")
        
        # Sort points by profile x-coordinate to determine adjacency
        profile_x_coords = strike_dip_data['profile_x']
        sorted_indices = np.argsort(profile_x_coords)
        
        # Analyze wedges between adjacent points (in profile plane order)
        arcpy.AddMessage("\nAnalyzing wedge geometry between adjacent points...")
        wedge_data = self.analyze_wedges(strike_dip_data, sorted_indices)
        
        # Output projected dip information
        arcpy.AddMessage("\nProjected dip angles in profile plane:")
        for i in range(len(strike_dip_data['labels'])):
            att = strike_dip_data['projected_attitudes'][i]
            # Calculate apparent dip angle in profile plane
            # The attitude vector (x,y) represents the trace, angle from horizontal
            dip_angle = np.degrees(np.arctan2(abs(att['y']), abs(att['x'])))
            # If y is negative, it's dipping down, if positive it's dipping up in the trace
            # But we want the angle from horizontal
            apparent_dip = 90 - dip_angle if abs(att['y']) > abs(att['x']) else dip_angle
            arcpy.AddMessage(f"  {strike_dip_data['labels'][i]}: {apparent_dip:.1f}°")
        
        # Output stratigraphic height vectors
        arcpy.AddMessage("\nStratigraphic height vectors in profile plane (x=right, y=up):")
        for i in range(len(strike_dip_data['labels'])):
            strat_vec = strike_dip_data['strat_height_vectors'][i]
            is_overturned = strike_dip_data['overturned'][i]
            status = "OVERTURNED" if is_overturned else "UPRIGHT"
            arcpy.AddMessage(f"  {strike_dip_data['labels'][i]} ({status}): x={strat_vec['x']:+.4f}, y={strat_vec['y']:+.4f}")
        
        # Output parallel bed information
        arcpy.AddMessage("\nParallel bed analysis:")
        has_parallel = False
        for wedge in wedge_data:
            if wedge['type'] in ['parallel', 'invalid_parallel']:
                has_parallel = True
                idx1 = wedge['point1_idx']
                idx2 = wedge['point2_idx']
                label1 = strike_dip_data['labels'][idx1]
                label2 = strike_dip_data['labels'][idx2]
                strat_vec1 = strike_dip_data['strat_height_vectors'][idx1]
                strat_vec2 = strike_dip_data['strat_height_vectors'][idx2]
                
                arcpy.AddMessage(f"  {label1} and {label2}: PARALLEL beds")
                arcpy.AddMessage(f"    {label1} strat height vector: x={strat_vec1['x']:+.4f}, y={strat_vec1['y']:+.4f}")
                arcpy.AddMessage(f"    {label2} strat height vector: x={strat_vec2['x']:+.4f}, y={strat_vec2['y']:+.4f}")
                
                if wedge['type'] == 'parallel':
                    arcpy.AddMessage(f"    -> VALID: Consistent stratigraphic directions")
                else:
                    arcpy.AddMessage(f"    -> INVALID: Opposite stratigraphic directions")
        
        if not has_parallel:
            arcpy.AddMessage("  No parallel beds found")
        
        arcpy.AddMessage(f"\nAnalyzed {len(wedge_data)} wedges")
        
        # Create downplunge projection
        arcpy.AddMessage("\nCreating downplunge projection...")
        
        # Extract DEM cell centroids if requested for plotting OR if wedge raster output is requested
        dem_cell_data = None
        if plot_dem_cells or wedge_raster_output:
            if not plot_dem_cells:
                arcpy.AddMessage("  Extracting DEM cell centroids for wedge assignment raster...")
            else:
                arcpy.AddMessage("  Extracting DEM cell centroids...")
            dem_cell_data = self.extract_dem_cell_centroids(dem_raster, study_area)
            arcpy.AddMessage(f"  Extracted {len(dem_cell_data['x'])} DEM cells")
        
        self.create_downplunge_view(strike_dip_data, wedge_data, dem_cell_data, fold_trend, fold_plunge, 
                                     output_pdf, dem_raster, study_area, wedge_raster_output, 
                                     strike_dip_fc, output_strike_dip_fc, label_field, plot_dem_cells,
                                     strat_height_raster_output)
        arcpy.AddMessage(f"  Downplunge view saved to: {output_pdf}")
        
        arcpy.AddMessage("\n" + "=" * 60)
        arcpy.AddMessage("Processing complete!")
        arcpy.AddMessage("=" * 60)
        
        return

    def check_utm_projection(self, dataset, dataset_name):
        """Check if dataset is in UTM projection"""
        desc = arcpy.Describe(dataset)
        sr = desc.spatialReference
        
        if not sr.type == "Projected":
            raise ValueError(f"{dataset_name} must be in a projected coordinate system (UTM required)")
        
        if "UTM" not in sr.name.upper():
            arcpy.AddWarning(f"WARNING: {dataset_name} may not be in UTM (found: {sr.name})")
            arcpy.AddWarning("  The tool expects UTM coordinates for accurate calculations")
        else:
            arcpy.AddMessage(f"  {dataset_name}: {sr.name} ✓")

    def read_strike_dip_data(self, fc, strike_field, dip_field, overturned_field, label_field, strat_height_field, dem_raster):
        """Read strike, dip, and location data from feature class"""
        
        # Check if feature class is Z-enabled
        desc = arcpy.Describe(fc)
        has_z = desc.hasZ
        
        fields = ["SHAPE@X", "SHAPE@Y", strike_field, dip_field, label_field]
        if has_z:
            fields.insert(2, "SHAPE@Z")
        if overturned_field:
            fields.append(overturned_field)
        if strat_height_field:
            fields.append(strat_height_field)
        
        data = {
            'x': [],
            'y': [],
            'z': [],
            'strike': [],
            'dip': [],
            'overturned': [],
            'labels': [],
            'strat_height': []
        }
        
        # Calculate field indices based on what's present
        current_idx = 0
        shape_x_idx = current_idx
        current_idx += 1
        shape_y_idx = current_idx
        current_idx += 1
        
        if has_z:
            z_index = current_idx
            current_idx += 1
        else:
            z_index = None
        
        strike_index = current_idx
        current_idx += 1
        dip_index = current_idx
        current_idx += 1
        label_index = current_idx
        current_idx += 1
        
        if overturned_field:
            overturned_index = current_idx
            current_idx += 1
        else:
            overturned_index = None
        
        if strat_height_field:
            strat_height_index = current_idx
            current_idx += 1
        else:
            strat_height_index = None
        
        with arcpy.da.SearchCursor(fc, fields) as cursor:
            for row in cursor:
                data['x'].append(row[shape_x_idx])
                data['y'].append(row[shape_y_idx])
                
                # Get Z from feature or set to None for later extraction from DEM
                if has_z and row[z_index] is not None:
                    data['z'].append(row[z_index])
                else:
                    data['z'].append(None)
                
                data['strike'].append(row[strike_index])
                data['dip'].append(row[dip_index])
                data['labels'].append(str(row[label_index]))
                
                # Check overturned flag
                if overturned_field:
                    data['overturned'].append(1 if row[overturned_index] == 1 else 0)
                else:
                    data['overturned'].append(0)
                
                # Get stratigraphic height if field provided
                if strat_height_field:
                    data['strat_height'].append(row[strat_height_index])
                else:
                    data['strat_height'].append(None)
        
        # Convert to numpy arrays
        for key in data:
            data[key] = np.array(data[key])
        
        # Extract Z values from DEM where needed
        if None in data['z']:
            arcpy.AddMessage("  Extracting elevations from DEM...")
            data['z'] = self.extract_z_from_dem(data['x'], data['y'], data['z'], dem_raster)
        
        return data
    
    def extract_z_from_dem(self, x_coords, y_coords, z_coords, dem_raster):
        """Extract Z values from DEM for points"""
        
        # Create temporary point array
        points = []
        for i in range(len(x_coords)):
            if z_coords[i] is None:
                points.append((x_coords[i], y_coords[i]))
        
        # Use arcpy to extract values
        point_geoms = [arcpy.PointGeometry(arcpy.Point(x, y)) for x, y in points]
        
        # Extract Multi Values to Points would be ideal but requires feature class
        # Instead, use Get Cell Value for each point
        z_values = []
        z_idx = 0
        for i in range(len(x_coords)):
            if z_coords[i] is None:
                # Get cell value at this location
                result = arcpy.management.GetCellValue(dem_raster, f"{x_coords[i]} {y_coords[i]}")
                cell_value = result.getOutput(0)
                try:
                    z_values.append(float(cell_value))
                except:
                    arcpy.AddWarning(f"  Warning: Could not extract DEM value at point {i+1}, using 0.0")
                    z_values.append(0.0)
            else:
                z_values.append(z_coords[i])
        
        return np.array(z_values)

    def calculate_fold_axis(self, strike_dip_data):
        """Calculate best-fit cylindrical fold axis from bedding measurements"""
        
        strikes = strike_dip_data['strike']
        dips = strike_dip_data['dip']
        
        # Convert strikes and dips to poles (unit vectors normal to bedding)
        poles = []
        for i in range(len(strikes)):
            strike = strikes[i]
            dip = dips[i]
            
            # Calculate pole to bedding
            # Pole trend = strike - 90° (perpendicular to strike, in dip direction)
            # Pole plunge = 90° - dip
            pole_trend = strike - 90
            pole_plunge = 90 - dip
            
            # Normalize trend to 0-360
            pole_trend = pole_trend % 360
            
            # Convert to Cartesian coordinates
            pole_trend_rad = np.radians(pole_trend)
            pole_plunge_rad = np.radians(pole_plunge)
            
            # Unit vector: x=E, y=N, z=Up
            px = np.cos(pole_plunge_rad) * np.sin(pole_trend_rad)
            py = np.cos(pole_plunge_rad) * np.cos(pole_trend_rad)
            pz = np.sin(pole_plunge_rad)
            
            poles.append([px, py, pz])
        
        poles = np.array(poles)
        
        # The fold axis is the eigenvector with smallest eigenvalue of poles^T * poles
        # This minimizes the sum of squared angles between poles and a plane perpendicular to the axis
        covariance = np.dot(poles.T, poles)
        eigenvalues, eigenvectors = np.linalg.eig(covariance)
        
        # Find eigenvector with smallest eigenvalue
        min_idx = np.argmin(eigenvalues)
        fold_axis = eigenvectors[:, min_idx]
        
        # Convert back to trend and plunge
        x, y, z = fold_axis
        
        # Plunge (can be negative initially)
        plunge = np.degrees(np.arcsin(z))
        
        # Trend (handle all quadrants correctly)
        trend = np.degrees(np.arctan2(x, y))
        if trend < 0:
            trend += 360
        
        # Convention: trend should point in direction of descent (plunge direction)
        # If plunge is negative, flip to make it positive and adjust trend by 180°
        if plunge < 0:
            plunge = -plunge
            trend = (trend + 180) % 360
        
        return trend, plunge

    def calculate_projected_attitudes(self, strike_dip_data, fold_trend, fold_plunge):
        """Calculate the intersection of each bedding plane with the profile plane"""
        
        strikes = strike_dip_data['strike']
        dips = strike_dip_data['dip']
        overturned = strike_dip_data['overturned']
        
        # Convert fold axis to unit vector
        trend_rad = np.radians(fold_trend)
        plunge_rad = np.radians(fold_plunge)
        
        fold_axis = np.array([
            np.cos(plunge_rad) * np.sin(trend_rad),  # x (East)
            np.cos(plunge_rad) * np.cos(trend_rad),  # y (North)
            np.sin(plunge_rad)                        # z (Up)
        ])
        
        # Define profile plane basis vectors
        vertical = np.array([0, 0, 1])
        v_parallel = np.dot(vertical, fold_axis) * fold_axis
        profile_y = vertical - v_parallel
        profile_y = profile_y / np.linalg.norm(profile_y)
        
        profile_x = np.cross(profile_y, fold_axis)
        profile_x = profile_x / np.linalg.norm(profile_x)
        
        # Project measurement points onto profile plane
        points_3d = np.column_stack([
            strike_dip_data['x'],
            strike_dip_data['y'],
            strike_dip_data['z']
        ])
        
        # Center the points
        centroid = np.mean(points_3d, axis=0)
        points_centered = points_3d - centroid
        
        # Project onto profile plane
        profile_x_coords = np.dot(points_centered, profile_x)
        profile_y_coords = np.dot(points_centered, profile_y)
        
        # Calculate pole to each bedding plane and its intersection with profile plane
        projected_attitudes = []
        
        for i in range(len(strikes)):
            strike = strikes[i]
            dip = dips[i]
            is_overturned = overturned[i]
            
            # Convert to radians
            strike_rad = np.radians(strike)
            dip_rad = np.radians(dip)
            
            # Calculate pole to bedding (right-hand rule)
            pole_trend = strike + 90
            pole_plunge = 90 - dip
            
            # If overturned, flip the pole
            if is_overturned:
                pole_trend += 180
                pole_plunge = -pole_plunge
            
            pole_trend = pole_trend % 360
            pole_trend_rad = np.radians(pole_trend)
            pole_plunge_rad = np.radians(pole_plunge)
            
            # Pole unit vector
            pole = np.array([
                np.cos(pole_plunge_rad) * np.sin(pole_trend_rad),
                np.cos(pole_plunge_rad) * np.cos(pole_trend_rad),
                np.sin(pole_plunge_rad)
            ])
            
            # The bedding plane has normal vector = pole
            # The profile plane has normal vector = fold_axis
            # The intersection line is perpendicular to both normals
            intersection_3d = np.cross(pole, fold_axis)
            
            # Normalize
            if np.linalg.norm(intersection_3d) > 1e-10:
                intersection_3d = intersection_3d / np.linalg.norm(intersection_3d)
            else:
                # Bedding plane is parallel to profile plane - use profile_x as default
                intersection_3d = profile_x.copy()
            
            # Project intersection line onto profile plane coordinates
            intersection_x = np.dot(intersection_3d, profile_x)
            intersection_y = np.dot(intersection_3d, profile_y)
            
            # Store as unit vector in profile plane
            norm = np.sqrt(intersection_x**2 + intersection_y**2)
            if norm > 1e-10:
                intersection_x /= norm
                intersection_y /= norm
            
            projected_attitudes.append({
                'x': intersection_x,
                'y': intersection_y,
                'pole_3d': pole,
                'intersection_3d': intersection_3d
            })
        
        strike_dip_data['projected_attitudes'] = projected_attitudes
        
        # Calculate stratigraphic height vectors in profile plane for each point
        # This is the perpendicular to bedding, oriented based on overturned flag
        strat_height_vectors = []
        
        # Also store profile coordinates for later use
        strike_dip_data['profile_x'] = profile_x_coords
        strike_dip_data['profile_y'] = profile_y_coords
        
        for i in range(len(strikes)):
            att = projected_attitudes[i]
            is_overturned = overturned[i]
            
            # Calculate perpendicular to bedding trace in profile plane
            perp_x = -att['y']
            perp_y = att['x']
            norm = np.sqrt(perp_x**2 + perp_y**2)
            if norm > 1e-10:
                perp_x /= norm
                perp_y /= norm
            
            # Apply overturned flag:
            # - Upright beds: perpendicular should point UP (positive y)
            # - Overturned beds: perpendicular should point DOWN (negative y)
            if is_overturned:
                # For overturned beds, ensure y is negative
                if perp_y > 0:
                    perp_x = -perp_x
                    perp_y = -perp_y
            else:
                # For upright beds, ensure y is positive
                if perp_y < 0:
                    perp_x = -perp_x
                    perp_y = -perp_y
            
            strat_height_vectors.append({'x': perp_x, 'y': perp_y})
        
        strike_dip_data['strat_height_vectors'] = strat_height_vectors
        
        return strike_dip_data
    
    def analyze_wedges(self, strike_dip_data, sorted_indices):
        """Analyze wedge geometry between adjacent measurement pairs (sorted by profile x)"""
        
        # Tolerance for determining if dips are parallel (in degrees)
        PARALLEL_TOLERANCE = 0.01
        
        wedge_data = []
        
        n_points = len(sorted_indices)
        
        # Create wedges between adjacent points in sorted order
        for i in range(n_points - 1):
            # Get original indices from sorted order
            idx1 = sorted_indices[i]
            idx2 = sorted_indices[i + 1]
            
            wedge = {
                'point1_idx': idx1,
                'point2_idx': idx2,
                'type': None,  # 'parallel' or 'intersecting'
                'intersection_point': None,
                'cylinder_axis_location': None
            }
            
            # Get projected attitudes for both points
            att1 = strike_dip_data['projected_attitudes'][idx1]
            att2 = strike_dip_data['projected_attitudes'][idx2]
            
            # Check if attitudes are parallel
            # Calculate angle between the two attitude vectors
            dot_product = att1['x'] * att2['x'] + att1['y'] * att2['y']
            # Clamp to [-1, 1] to avoid numerical issues with arccos
            dot_product = np.clip(dot_product, -1.0, 1.0)
            angle_diff = np.degrees(np.arccos(np.abs(dot_product)))
            
            if angle_diff < PARALLEL_TOLERANCE:
                # Check if stratigraphic height vectors point in same geometric direction
                strat_vec1 = strike_dip_data['strat_height_vectors'][idx1]
                strat_vec2 = strike_dip_data['strat_height_vectors'][idx2]
                
                # If y components have opposite signs, they're geometrically opposite
                if strat_vec1['y'] * strat_vec2['y'] >= 0:
                    # Same geometric direction - valid parallel case
                    wedge['type'] = 'parallel'
                else:
                    # Opposite directions - treat as invalid
                    wedge['type'] = 'invalid_parallel'
            
            wedge_data.append(wedge)
        
        return wedge_data
    
    def extract_dem_cell_centroids(self, dem_raster, study_area):
        """Extract centroids and elevation values for each DEM cell within study area"""
        
        # Get DEM properties
        dem_desc = arcpy.Describe(dem_raster)
        cell_size_x = dem_desc.meanCellWidth
        cell_size_y = dem_desc.meanCellHeight
        extent = dem_desc.extent
        
        # Clip DEM to study area if provided
        if study_area:
            arcpy.AddMessage("    Clipping DEM to study area...")
            clipped_dem = arcpy.sa.ExtractByMask(dem_raster, study_area)
        else:
            clipped_dem = dem_raster
        
        # Convert raster to numpy array
        arcpy.AddMessage("    Converting DEM to point array...")
        
        # Get NoData value from the original or clipped raster
        try:
            nodata_value = arcpy.Describe(clipped_dem).noDataValue
        except:
            # If clipped_dem is a Raster object, get from original
            nodata_value = arcpy.Describe(dem_raster).noDataValue
        
        if nodata_value is None:
            nodata_value = -9999
        
        # Convert to numpy array
        dem_array = arcpy.RasterToNumPyArray(clipped_dem, nodata_to_value=nodata_value)
        
        # Get extent for coordinate calculation
        if study_area:
            # For Raster objects from ExtractByMask, we need to get extent differently
            try:
                clipped_extent = arcpy.Describe(clipped_dem).extent
            except:
                # Save temporarily to get extent
                import tempfile
                import os
                temp_raster = os.path.join(tempfile.gettempdir(), "temp_clipped.tif")
                clipped_dem.save(temp_raster)
                clipped_extent = arcpy.Describe(temp_raster).extent
                arcpy.Delete_management(temp_raster)
            
            ll_corner = clipped_extent.lowerLeft
            extent = clipped_extent
        else:
            ll_corner = extent.lowerLeft
        
        # Generate cell centroids
        cell_data = {
            'x': [],
            'y': [],
            'z': []
        }
        
        nrows, ncols = dem_array.shape
        
        for row in range(nrows):
            for col in range(ncols):
                elev = dem_array[row, col]
                
                # Skip NoData cells
                if elev == nodata_value:
                    continue
                
                # Calculate centroid coordinates
                # Note: array[0,0] is upper-left, so need to flip row indexing
                x = ll_corner.X + (col + 0.5) * cell_size_x
                y = ll_corner.Y + (nrows - row - 0.5) * cell_size_y
                
                cell_data['x'].append(x)
                cell_data['y'].append(y)
                cell_data['z'].append(elev)
        
        # Convert to numpy arrays
        for key in cell_data:
            cell_data[key] = np.array(cell_data[key])
        
        return cell_data
    
    def assign_dem_cells_to_wedges(self, dem_profile_x, dem_profile_y, strike_dip_data, wedge_data, sorted_indices):
        """Assign each DEM cell to a wedge or rectangle"""
        
        arcpy.AddMessage("  Assigning DEM cells to wedges/rectangles...")
        
        # Initialize wedge assignments (-1 = not in any wedge, -2 = in multiple wedges)
        wedge_assignments = np.full(len(dem_profile_x), -1, dtype=int)
        
        # Set up progress tracking
        n_wedges = len([w for w in wedge_data if w['type'] in ['parallel', 'intersecting']])
        arcpy.SetProgressor("step", "Assigning DEM cells to wedges...", 0, n_wedges, 1)
        
        wedge_counter = 0
        
        # Process each wedge in order
        for wedge_idx, wedge in enumerate(wedge_data):
            if wedge['type'] not in ['parallel', 'intersecting']:
                continue
            
            wedge_counter += 1
            percent = int(100 * wedge_counter / n_wedges)
            arcpy.SetProgressorLabel(f"Assigning DEM cells to wedges... {percent}%")
            arcpy.SetProgressorPosition(wedge_counter)
            
            left_idx = wedge['point1_idx']
            right_idx = wedge['point2_idx']
            
            # Get positions in profile plane
            left_x = strike_dip_data['profile_x'][left_idx]
            left_y = strike_dip_data['profile_y'][left_idx]
            right_x = strike_dip_data['profile_x'][right_idx]
            right_y = strike_dip_data['profile_y'][right_idx]
            
            # Get strat height vectors
            strat_vec_left = strike_dip_data['strat_height_vectors'][left_idx]
            strat_vec_right = strike_dip_data['strat_height_vectors'][right_idx]
            
            if wedge['type'] == 'parallel':
                # RECTANGLE CASE
                # Translate so left point is at origin
                dem_translated_x = dem_profile_x - left_x
                dem_translated_y = dem_profile_y - left_y
                right_translated_x = right_x - left_x
                right_translated_y = right_y - left_y
                
                # Rotate so left strat height vector points along +y axis
                # angle is angle of strat vector from +x axis; rotate by (pi/2 - angle) to bring onto +y
                angle = np.arctan2(strat_vec_left['y'], strat_vec_left['x'])
                rot_angle = np.pi/2 - angle
                cos_a = np.cos(rot_angle)
                sin_a = np.sin(rot_angle)
                
                # Rotate DEM points
                dem_rot_x = cos_a * dem_translated_x - sin_a * dem_translated_y
                dem_rot_y = sin_a * dem_translated_x + cos_a * dem_translated_y
                
                # Rotate right point
                right_rot_x = cos_a * right_translated_x - sin_a * right_translated_y
                right_rot_y = sin_a * right_translated_x + cos_a * right_translated_y
                
                # Check if DEM point is in rectangle
                # Point is inside if between 0 and right point in both x and y
                if right_rot_x >= 0:
                    in_x = (dem_rot_x >= 0) & (dem_rot_x <= right_rot_x)
                else:
                    in_x = (dem_rot_x <= 0) & (dem_rot_x >= right_rot_x)
                
                if right_rot_y >= 0:
                    in_y = (dem_rot_y >= 0) & (dem_rot_y <= right_rot_y)
                else:
                    in_y = (dem_rot_y <= 0) & (dem_rot_y >= right_rot_y)
                
                in_rectangle = in_x & in_y
                
                # Mark points: if already assigned, set to -2 (multiple wedges)
                already_assigned = (wedge_assignments >= 0) & in_rectangle
                wedge_assignments[already_assigned] = -2
                # Assign unassigned points
                newly_assigned = (wedge_assignments == -1) & in_rectangle
                wedge_assignments[newly_assigned] = wedge_idx
                
            elif wedge['type'] == 'intersecting':
                # WEDGE CASE
                # Translate so intersection point is at origin
                intersection = wedge['intersection_point']
                dem_translated_x = dem_profile_x - intersection[0]
                dem_translated_y = dem_profile_y - intersection[1]
                
                # Get measurement point positions
                left_x = strike_dip_data['profile_x'][left_idx]
                left_y = strike_dip_data['profile_y'][left_idx]
                right_x = strike_dip_data['profile_x'][right_idx]
                right_y = strike_dip_data['profile_y'][right_idx]
                
                # Calculate vectors FROM intersection TO measurement points
                # These define the wedge boundaries
                boundary1 = np.array([left_x - intersection[0], left_y - intersection[1]])
                boundary2 = np.array([right_x - intersection[0], right_y - intersection[1]])
                
                # Normalize
                boundary1 = boundary1 / np.linalg.norm(boundary1)
                boundary2 = boundary2 / np.linalg.norm(boundary2)
                
                # Calculate dot product between boundaries
                dot_between = np.dot(boundary1, boundary2)
                
                # For each DEM point, check if inside wedge
                for i in range(len(dem_profile_x)):
                    # Vector to DEM point from intersection
                    dem_vec = np.array([dem_translated_x[i], dem_translated_y[i]])
                    dem_dist = np.linalg.norm(dem_vec)
                    
                    if dem_dist < 1e-10:
                        # Point is at the intersection - skip
                        continue
                    
                    # Normalize
                    dem_vec = dem_vec / dem_dist
                    
                    # Calculate dot products
                    dot1 = np.dot(boundary1, dem_vec)
                    dot2 = np.dot(boundary2, dem_vec)
                    
                    # Point is inside wedge if both dots are greater than dot_between
                    if dot1 > dot_between and dot2 > dot_between:
                        if wedge_assignments[i] >= 0:
                            # Already assigned to another wedge
                            wedge_assignments[i] = -2
                        elif wedge_assignments[i] == -1:
                            # First assignment
                            wedge_assignments[i] = wedge_idx
        
        n_assigned = np.sum(wedge_assignments >= 0)
        n_multiple = np.sum(wedge_assignments == -2)
        
        arcpy.ResetProgressor()
        
        arcpy.AddMessage(f"  Assigned {n_assigned} DEM cells to wedges/rectangles")
        if n_multiple > 0:
            arcpy.AddMessage(f"  Warning: {n_multiple} DEM cells belong to multiple wedges (marked as ambiguous)")
        
        return wedge_assignments
    
    def calculate_dem_stratigraphic_heights(self, dem_profile_x, dem_profile_y, dem_wedge_assignments, 
                                           strike_dip_data, wedge_data):
        """Calculate stratigraphic heights for DEM cells based on their wedge assignments"""
        
        arcpy.AddMessage("  Calculating stratigraphic heights for DEM cells...")
        
        # Initialize heights array (NaN for unassigned)
        dem_strat_heights = np.full(len(dem_profile_x), np.nan)
        
        # Get calculated heights for measurement points
        if 'calculated_strat_height' not in strike_dip_data:
            arcpy.AddMessage("  Warning: No calculated stratigraphic heights available for measurements")
            return dem_strat_heights
        
        point_heights = strike_dip_data['calculated_strat_height']
        
        # Set up progress tracking
        n_wedges = len([w for w in wedge_data if w['type'] in ['parallel', 'intersecting']])
        arcpy.SetProgressor("step", "Calculating DEM stratigraphic heights...", 0, n_wedges, 1)
        
        wedge_counter = 0
        
        # Process each wedge
        for wedge_idx, wedge in enumerate(wedge_data):
            if wedge['type'] not in ['parallel', 'intersecting']:
                continue
            
            wedge_counter += 1
            percent = int(100 * wedge_counter / n_wedges)
            arcpy.SetProgressorLabel(f"Calculating DEM stratigraphic heights... {percent}%")
            arcpy.SetProgressorPosition(wedge_counter)
            
            # Get indices of DEM cells in this wedge
            mask = (dem_wedge_assignments == wedge_idx)
            n_cells = np.sum(mask)
            if n_cells == 0:
                continue
            
            left_idx = wedge['point1_idx']
            right_idx = wedge['point2_idx']
            
            # Get left point height (reference)
            left_height = point_heights[left_idx]
            if left_height is None:
                arcpy.AddWarning(f"  Warning: No height available for left point in wedge {wedge_idx}")
                continue
            
            # Get positions in profile plane
            left_x = strike_dip_data['profile_x'][left_idx]
            left_y = strike_dip_data['profile_y'][left_idx]
            
            if wedge['type'] == 'parallel':
                # RECTANGLE CASE
                # Rotate coordinate system so left strat height vector points along +y
                strat_vec_left = strike_dip_data['strat_height_vectors'][left_idx]
                
                # Get DEM points in this wedge
                dem_x = dem_profile_x[mask]
                dem_y = dem_profile_y[mask]
                
                # Translate so left point is at origin
                dem_translated_x = dem_x - left_x
                dem_translated_y = dem_y - left_y
                
                # Rotate so left strat height vector points along +y axis
                # angle is angle of strat vector from +x axis; rotate by (pi/2 - angle) to bring onto +y
                angle = np.arctan2(strat_vec_left['y'], strat_vec_left['x'])
                rot_angle = np.pi/2 - angle
                cos_a = np.cos(rot_angle)
                sin_a = np.sin(rot_angle)
                
                dem_rot_x = cos_a * dem_translated_x - sin_a * dem_translated_y
                dem_rot_y = sin_a * dem_translated_x + cos_a * dem_translated_y
                
                # Stratigraphic height = y coordinate in rotated system + left point height
                heights = dem_rot_y + left_height
                dem_strat_heights[mask] = heights
                
            elif wedge['type'] == 'intersecting':
                # WEDGE CASE
                intersection = wedge['intersection_point']
                intersection_type = wedge.get('intersection_type', 'unknown')
                
                # Get DEM points in this wedge
                dem_x = dem_profile_x[mask]
                dem_y = dem_profile_y[mask]
                
                # Distance from each DEM point to intersection
                dx_dem = dem_x - intersection[0]
                dy_dem = dem_y - intersection[1]
                dist_from_intersection = np.sqrt(dx_dem**2 + dy_dem**2)
                
                # Distance from left point to intersection
                dx_left = left_x - intersection[0]
                dy_left = left_y - intersection[1]
                dist_left_to_intersection = np.sqrt(dx_left**2 + dy_left**2)
                
                if intersection_type == 'up-up':
                    # strat_height = -distance_from_intersection + constant
                    # constant = left_height + dist_left_to_intersection
                    constant = left_height + dist_left_to_intersection
                    heights = -dist_from_intersection + constant
                    
                elif intersection_type == 'down-down':
                    # strat_height = +distance_from_intersection + constant
                    # constant = left_height - dist_left_to_intersection
                    constant = left_height - dist_left_to_intersection
                    heights = dist_from_intersection + constant
                    
                else:
                    arcpy.AddWarning(f"  Warning: Unknown intersection type for wedge {wedge_idx}")
                    continue
                
                dem_strat_heights[mask] = heights
        
        n_calculated = np.sum(~np.isnan(dem_strat_heights))
        
        arcpy.ResetProgressor()
        
        arcpy.AddMessage(f"  Calculated stratigraphic heights for {n_calculated} DEM cells")
        
        return dem_strat_heights
    
    def create_wedge_assignment_raster(self, dem_cell_data, wedge_assignments, dem_raster, study_area, output_path):
        """Create a raster showing which wedge each DEM cell belongs to"""
        
        # Get DEM properties for creating output raster
        dem_desc = arcpy.Describe(dem_raster)
        cell_size_x = dem_desc.meanCellWidth
        cell_size_y = dem_desc.meanCellHeight
        
        # Clip DEM to study area to match dimensions
        if study_area:
            clipped_dem = arcpy.sa.ExtractByMask(dem_raster, study_area)
        else:
            clipped_dem = dem_raster
        
        # Get extent and nodata
        try:
            clipped_extent = arcpy.Describe(clipped_dem).extent
            nodata_value = arcpy.Describe(clipped_dem).noDataValue
        except:
            # If clipped_dem is a Raster object, save temporarily
            import tempfile
            import os
            temp_raster = os.path.join(tempfile.gettempdir(), "temp_wedge.tif")
            clipped_dem.save(temp_raster)
            clipped_extent = arcpy.Describe(temp_raster).extent
            nodata_value = arcpy.Describe(temp_raster).noDataValue
            arcpy.Delete_management(temp_raster)
        
        if nodata_value is None:
            nodata_value = -9999
        
        # Convert DEM to array to get dimensions
        dem_array = arcpy.RasterToNumPyArray(clipped_dem, nodata_to_value=nodata_value)
        nrows, ncols = dem_array.shape
        
        # Create output array initialized with nodata
        wedge_array = np.full((nrows, ncols), -9999, dtype=np.int16)
        
        # Fill in wedge assignments
        ll_corner = clipped_extent.lowerLeft
        
        for i in range(len(dem_cell_data['x'])):
            x = dem_cell_data['x'][i]
            y = dem_cell_data['y'][i]
            wedge_id = wedge_assignments[i]
            
            # Calculate row and column
            col = int((x - ll_corner.X) / cell_size_x)
            row = int((ll_corner.Y + nrows * cell_size_y - y) / cell_size_y)
            
            # Check bounds
            if 0 <= row < nrows and 0 <= col < ncols:
                wedge_array[row, col] = wedge_id
        
        # Convert array to raster
        wedge_raster = arcpy.NumPyArrayToRaster(
            wedge_array,
            lower_left_corner=ll_corner,
            x_cell_size=cell_size_x,
            y_cell_size=cell_size_y,
            value_to_nodata=-9999
        )
        
        # Set spatial reference
        arcpy.management.DefineProjection(wedge_raster, dem_desc.spatialReference)
        
        # Save the raster
        wedge_raster.save(output_path)
    
    def calculate_stratigraphic_heights_in_profile(self, strike_dip_data, wedge_data, profile_x_coords, profile_y_coords, sorted_indices):
        """Calculate stratigraphic heights for measurements using the Busk method"""
        
        # Validate that exactly one measurement has a known stratigraphic height
        strat_heights = strike_dip_data['strat_height']
        known_heights = [h for h in strat_heights if h is not None]
        
        if len(known_heights) == 0:
            arcpy.AddWarning("No stratigraphic height specified - cannot calculate heights")
            return
        elif len(known_heights) > 1:
            arcpy.AddWarning(f"Multiple measurements have stratigraphic heights specified ({len(known_heights)} found)")
            arcpy.AddWarning("Only one measurement should have a known stratigraphic height - cannot calculate heights")
            return
        
        # Find the measurement with known height
        known_idx = None
        known_height = None
        for i, h in enumerate(strat_heights):
            if h is not None:
                known_idx = i
                known_height = h
                break
        
        arcpy.AddMessage(f"  Using {strike_dip_data['labels'][known_idx]} as reference with height = {known_height:.2f} m")
        
        # Initialize all heights as None, will calculate starting from leftmost = 0
        calculated_heights = [None] * len(strike_dip_data['x'])
        
        # Start propagation from leftmost point (sorted_indices[0]) with height = 0
        leftmost_idx = sorted_indices[0]
        calculated_heights[leftmost_idx] = 0.0
        arcpy.AddMessage(f"  Starting propagation from leftmost point {strike_dip_data['labels'][leftmost_idx]} with initial height = 0.0")
        
        # Propagate through wedges from left to right
        for i in range(len(sorted_indices) - 1):
            left_idx = sorted_indices[i]
            right_idx = sorted_indices[i + 1]
            
            # Find the wedge between these two points
            wedge = None
            for w in wedge_data:
                if (w['point1_idx'] == left_idx and w['point2_idx'] == right_idx) or \
                   (w['point1_idx'] == right_idx and w['point2_idx'] == left_idx):
                    wedge = w
                    break
            
            if wedge is None or wedge['type'] not in ['parallel', 'intersecting']:
                arcpy.AddWarning(f"  No valid wedge found between {strike_dip_data['labels'][left_idx]} and {strike_dip_data['labels'][right_idx]}")
                continue
            
            left_height = calculated_heights[left_idx]
            left_x = profile_x_coords[left_idx]
            left_y = profile_y_coords[left_idx]
            right_x = profile_x_coords[right_idx]
            right_y = profile_y_coords[right_idx]
            
            if wedge['type'] == 'intersecting':
                # Get intersection point and type
                ix, iy = wedge['intersection_point']
                intersection_type = wedge.get('intersection_type', None)
                
                # Calculate distances from intersection
                dist_left = np.sqrt((left_x - ix)**2 + (left_y - iy)**2)
                dist_right = np.sqrt((right_x - ix)**2 + (right_y - iy)**2)
                
                arcpy.AddMessage(f"    Intersection at ({ix:.2f}, {iy:.2f})")
                arcpy.AddMessage(f"    Left point ({left_x:.2f}, {left_y:.2f}), dist={dist_left:.2f}")
                arcpy.AddMessage(f"    Right point ({right_x:.2f}, {right_y:.2f}), dist={dist_right:.2f}")
                
                if intersection_type == 'up-up':
                    # Up-up wedge: strat up points away from intersection
                    # At intersection = youngest rocks
                    # Moving away from intersection = moving down-section (older)
                    # height = -distance + constant
                    constant = left_height + dist_left
                    right_height = -dist_right + constant
                    arcpy.AddMessage(f"    UP-UP wedge: {strike_dip_data['labels'][left_idx]} ({left_height:.2f}) -> {strike_dip_data['labels'][right_idx]} ({right_height:.2f})")
                    
                elif intersection_type == 'down-down':
                    # Down-down wedge: strat down points toward intersection
                    # At intersection = oldest rocks
                    # Moving away from intersection = moving up-section (younger)
                    # height = +distance + constant
                    constant = left_height - dist_left
                    right_height = dist_right + constant
                    arcpy.AddMessage(f"    DOWN-DOWN wedge: {strike_dip_data['labels'][left_idx]} ({left_height:.2f}) -> {strike_dip_data['labels'][right_idx]} ({right_height:.2f})")
                    
                else:
                    arcpy.AddWarning(f"    Unknown intersection type for wedge between {strike_dip_data['labels'][left_idx]} and {strike_dip_data['labels'][right_idx]}")
                    continue
                
                calculated_heights[right_idx] = right_height
                
            elif wedge['type'] == 'parallel':
                # For parallel beds (rectangle)
                # Project the position of right point onto the strat height vector of left point
                strat_vec_left = strike_dip_data['strat_height_vectors'][left_idx]
                
                # Vector from left to right point
                dx = right_x - left_x
                dy = right_y - left_y
                
                # Project this displacement onto the strat height vector
                # projection = (displacement · strat_vec) = component along strat height direction
                projection = dx * strat_vec_left['x'] + dy * strat_vec_left['y']
                
                # The stratigraphic height change is this projection
                right_height = left_height + projection
                
                if projection > 0:
                    arcpy.AddMessage(f"    PARALLEL (up): {strike_dip_data['labels'][left_idx]} ({left_height:.2f}) -> {strike_dip_data['labels'][right_idx]} ({right_height:.2f}), offset={projection:.2f}")
                else:
                    arcpy.AddMessage(f"    PARALLEL (down): {strike_dip_data['labels'][left_idx]} ({left_height:.2f}) -> {strike_dip_data['labels'][right_idx]} ({right_height:.2f}), offset={projection:.2f}")
                
                calculated_heights[right_idx] = right_height
        
        # Calculate correction constant based on known height
        if calculated_heights[known_idx] is not None:
            correction = known_height - calculated_heights[known_idx]
            arcpy.AddMessage(f"\n  Applying correction constant: {correction:.2f} m")
            
            # Apply correction to all calculated heights
            for i in range(len(calculated_heights)):
                if calculated_heights[i] is not None:
                    calculated_heights[i] += correction
            
            # Output final stratigraphic heights
            arcpy.AddMessage(f"\n  Final stratigraphic heights:")
            for i in range(len(calculated_heights)):
                if calculated_heights[i] is not None:
                    arcpy.AddMessage(f"    {strike_dip_data['labels'][i]}: {calculated_heights[i]:.2f} m")
        else:
            arcpy.AddWarning(f"  Could not propagate to reference point {strike_dip_data['labels'][known_idx]}")
        
        # Store results
        strike_dip_data['calculated_strat_height'] = calculated_heights
    
    def create_output_strike_dip_fc(self, input_fc, strike_dip_data, output_fc, label_field):
        """Create output feature class with calculated stratigraphic heights added"""
        
        arcpy.AddMessage(f"  Creating output feature class: {output_fc}")
        
        # Get input feature class properties
        desc = arcpy.Describe(input_fc)
        
        # Parse output path
        import os
        output_path = os.path.dirname(output_fc)
        output_name = os.path.basename(output_fc)
        
        # Create new feature class with correct alias
        arcpy.AddMessage(f"  Creating feature class with alias: {output_name}")
        arcpy.management.CreateFeatureclass(
            out_path=output_path,
            out_name=output_name,
            geometry_type=desc.shapeType,
            has_m="ENABLED" if desc.hasM else "DISABLED",
            has_z="ENABLED" if desc.hasZ else "DISABLED",
            spatial_reference=desc.spatialReference,
            out_alias=output_name  # Set alias to match output name
        )
        
        # Copy fields from input (except OID and Shape)
        input_fields = arcpy.ListFields(input_fc)
        for field in input_fields:
            if field.type not in ['OID', 'Geometry']:
                arcpy.management.AddField(
                    output_fc,
                    field.name,
                    field.type,
                    field_precision=field.precision,
                    field_scale=field.scale,
                    field_length=field.length,
                    field_alias=field.aliasName,
                    field_is_nullable=field.isNullable
                )
        
        # Copy data from input to output
        arcpy.AddMessage(f"  Copying data from input feature class...")
        field_names = [f.name for f in input_fields if f.type not in ['OID', 'Geometry']]
        field_names.insert(0, 'SHAPE@')
        
        with arcpy.da.SearchCursor(input_fc, field_names) as search_cursor:
            with arcpy.da.InsertCursor(output_fc, field_names) as insert_cursor:
                for row in search_cursor:
                    insert_cursor.insertRow(row)
        
        arcpy.AddMessage(f"  Data copied successfully")
        
        # Add the strat_height_out field
        arcpy.AddMessage(f"  Adding 'strat_height_out' field...")
        arcpy.management.AddField(output_fc, 'strat_height_out', 'DOUBLE', 
                                 field_alias='Calculated Stratigraphic Height')
        
        # Update the field with calculated values
        calculated_heights = strike_dip_data.get('calculated_strat_height', [])
        labels = strike_dip_data['labels']
        
        arcpy.AddMessage(f"  Found {len(calculated_heights)} calculated heights")
        arcpy.AddMessage(f"  Found {len(labels)} labels")
        arcpy.AddMessage(f"  Sample labels: {labels[:3] if len(labels) > 0 else 'None'}")
        arcpy.AddMessage(f"  Sample heights: {[h for h in calculated_heights[:3] if h is not None]}")
        arcpy.AddMessage(f"  Using label field: {label_field}")
        
        # Create a mapping from label to calculated height
        label_to_height = {}
        for i, label in enumerate(labels):
            if i < len(calculated_heights) and calculated_heights[i] is not None:
                label_to_height[str(label)] = calculated_heights[i]
                if i < 3:  # Only show first 3 for debugging
                    arcpy.AddMessage(f"  Mapping: label '{label}' -> height {calculated_heights[i]:.2f}")
        
        arcpy.AddMessage(f"  Created mapping for {len(label_to_height)} features with heights")
        
        # Update features using the actual label field
        update_count = 0
        with arcpy.da.UpdateCursor(output_fc, [label_field, 'strat_height_out']) as cursor:
            for row in cursor:
                label_value = row[0]
                label_str = str(label_value)
                
                if label_str in label_to_height:
                    row[1] = label_to_height[label_str]
                    cursor.updateRow(row)
                    update_count += 1
                    if update_count <= 3:
                        arcpy.AddMessage(f"  Updated feature with label '{label_value}' -> height {label_to_height[label_str]:.2f}")
        
        arcpy.AddMessage(f"  Updated {update_count} features with calculated heights")
    
    def create_strat_height_raster(self, dem_cell_data, dem_strat_heights, dem_raster, study_area, output_path):
        """Create a raster showing stratigraphic height for each DEM cell"""
        
        # Get DEM properties for creating output raster
        dem_desc = arcpy.Describe(dem_raster)
        cell_size_x = dem_desc.meanCellWidth
        cell_size_y = dem_desc.meanCellHeight
        
        # Clip DEM to study area to match dimensions
        if study_area:
            clipped_dem = arcpy.sa.ExtractByMask(dem_raster, study_area)
        else:
            clipped_dem = dem_raster
        
        # Get extent and nodata
        try:
            clipped_extent = arcpy.Describe(clipped_dem).extent
            nodata_value = arcpy.Describe(clipped_dem).noDataValue
        except:
            # If clipped_dem is a Raster object, save temporarily
            import tempfile
            import os
            temp_raster = os.path.join(tempfile.gettempdir(), "temp_strat.tif")
            clipped_dem.save(temp_raster)
            clipped_extent = arcpy.Describe(temp_raster).extent
            nodata_value = arcpy.Describe(temp_raster).noDataValue
            arcpy.Delete_management(temp_raster)
        
        if nodata_value is None:
            nodata_value = -9999
        
        # Convert DEM to array to get dimensions
        dem_array = arcpy.RasterToNumPyArray(clipped_dem, nodata_to_value=nodata_value)
        nrows, ncols = dem_array.shape
        
        # Create output array initialized with nodata (-9999 for areas without height)
        strat_array = np.full((nrows, ncols), -9999.0, dtype=np.float32)
        
        # Fill in stratigraphic heights
        ll_corner = clipped_extent.lowerLeft
        
        for i in range(len(dem_cell_data['x'])):
            x = dem_cell_data['x'][i]
            y = dem_cell_data['y'][i]
            strat_height = dem_strat_heights[i]
            
            # Skip if no height calculated (NaN)
            if np.isnan(strat_height):
                continue
            
            # Calculate row and column
            col = int((x - ll_corner.X) / cell_size_x)
            row = int((ll_corner.Y + nrows * cell_size_y - y) / cell_size_y)
            
            # Check bounds
            if 0 <= row < nrows and 0 <= col < ncols:
                strat_array[row, col] = strat_height
        
        # Convert array to raster
        strat_raster = arcpy.NumPyArrayToRaster(
            strat_array,
            lower_left_corner=ll_corner,
            x_cell_size=cell_size_x,
            y_cell_size=cell_size_y,
            value_to_nodata=-9999
        )
        
        # Set spatial reference
        arcpy.management.DefineProjection(strat_raster, dem_desc.spatialReference)
        
        # Save the raster
        strat_raster.save(output_path)

    def create_downplunge_view(self, strike_dip_data, wedge_data, dem_cell_data, fold_trend, fold_plunge, 
                                output_pdf, dem_raster, study_area, wedge_raster_output, 
                                strike_dip_fc, output_strike_dip_fc, label_field, plot_dem_cells,
                                strat_height_raster_output):
        """Create downplunge projection and save as PDF"""
        
        # Convert fold axis to unit vector
        # Trend = azimuth direction of plunge (0-360°)
        # Plunge = angle below horizontal (0-90°)
        # The fold axis points in the direction of plunge (downward and in trend direction)
        trend_rad = np.radians(fold_trend)
        plunge_rad = np.radians(fold_plunge)
        
        fold_axis = np.array([
            np.cos(plunge_rad) * np.sin(trend_rad),  # x (East component)
            np.cos(plunge_rad) * np.cos(trend_rad),  # y (North component)
            -np.sin(plunge_rad)                       # z (Down, negative for plunge)
        ])
        
        # For down-plunge view:
        # - Fold axis points in the viewing direction (into the page/screen)
        # - Profile y-axis points upward in the view
        # - Profile x-axis points to the right in the view
        
        vertical = np.array([0, 0, 1])
        
        # Project vertical onto profile plane (remove component parallel to fold axis)
        # This gives us the "up" direction in the profile view
        v_parallel = np.dot(vertical, fold_axis) * fold_axis
        profile_y = vertical - v_parallel
        profile_y = profile_y / np.linalg.norm(profile_y)
        
        # Profile x-axis: cross product to get "right" direction
        # Use profile_y × fold_axis to get right-hand direction
        profile_x = np.cross(profile_y, fold_axis)
        profile_x = profile_x / np.linalg.norm(profile_x)
        
        # Project measurement points onto profile plane
        points_3d = np.column_stack([
            strike_dip_data['x'],
            strike_dip_data['y'],
            strike_dip_data['z']
        ])
        
        # Center the points
        centroid = np.mean(points_3d, axis=0)
        points_centered = points_3d - centroid
        
        # Project onto profile plane
        profile_x_coords = np.dot(points_centered, profile_x)
        profile_y_coords = np.dot(points_centered, profile_y)
        
        # Update the stored profile coordinates (they should already match from calculate_projected_attitudes)
        strike_dip_data['profile_x'] = profile_x_coords
        strike_dip_data['profile_y'] = profile_y_coords
        
        # Sort points by profile x-coordinate (should already be done, but ensure it's here)
        sorted_indices = np.argsort(profile_x_coords)
        
        # Project DEM cells if provided
        dem_profile_x = None
        dem_profile_y = None
        dem_wedge_assignments = None
        if dem_cell_data is not None:
            arcpy.AddMessage("  Projecting DEM cells onto profile plane...")
            dem_points_3d = np.column_stack([
                dem_cell_data['x'],
                dem_cell_data['y'],
                dem_cell_data['z']
            ])
            dem_points_centered = dem_points_3d - centroid
            dem_profile_x = np.dot(dem_points_centered, profile_x)
            dem_profile_y = np.dot(dem_points_centered, profile_y)
        
        # Calculate wedge intersections now that we have profile coordinates
        intersection_points = []
        for wedge in wedge_data:
            if wedge['type'] == 'parallel':
                continue
            
            if wedge['type'] == 'invalid_parallel':
                # Parallel but opposite stratigraphic directions
                idx1 = wedge['point1_idx']
                idx2 = wedge['point2_idx']
                label1 = strike_dip_data['labels'][idx1]
                label2 = strike_dip_data['labels'][idx2]
                arcpy.AddWarning(f"No intersection between stratigraphic up or stratigraphic down for measurements {label1} and {label2}")
                continue
            
            # Get indices in original order
            idx1 = wedge['point1_idx']
            idx2 = wedge['point2_idx']
            
            label1 = strike_dip_data['labels'][idx1]
            label2 = strike_dip_data['labels'][idx2]
            
            # Find where these indices are in sorted order
            sorted_pos1 = np.where(sorted_indices == idx1)[0][0]
            sorted_pos2 = np.where(sorted_indices == idx2)[0][0]
            
            # Only process if they're adjacent in sorted order
            if abs(sorted_pos1 - sorted_pos2) != 1:
                continue
            
            # Get profile coordinates
            p1_x = profile_x_coords[idx1]
            p1_y = profile_y_coords[idx1]
            p2_x = profile_x_coords[idx2]
            p2_y = profile_y_coords[idx2]
            
            # Get the stratigraphic height vectors from the stored values
            strat_vec1 = strike_dip_data['strat_height_vectors'][idx1]
            strat_vec2 = strike_dip_data['strat_height_vectors'][idx2]
            
            # Debug output
            arcpy.AddMessage(f"\nChecking intersection for {label1} and {label2}:")
            arcpy.AddMessage(f"  {label1} strat height vector: x={strat_vec1['x']:+.4f}, y={strat_vec1['y']:+.4f}")
            arcpy.AddMessage(f"  {label2} strat height vector: x={strat_vec2['x']:+.4f}, y={strat_vec2['y']:+.4f}")
            
            # Find intersection of stratigraphic height vectors
            intersection, t1, t2 = self.find_line_intersection_with_params(
                p1_x, p1_y, strat_vec1['x'], strat_vec1['y'],
                p2_x, p2_y, strat_vec2['x'], strat_vec2['y']
            )
            
            if intersection is None:
                arcpy.AddMessage(f"  -> No intersection (parallel vectors)")
                arcpy.AddWarning(f"No intersection between stratigraphic height vectors for measurements {label1} and {label2}")
                continue
            
            arcpy.AddMessage(f"  -> Intersection at ({intersection[0]:.2f}, {intersection[1]:.2f})")
            arcpy.AddMessage(f"  -> Parameters: t1={t1:.2f}, t2={t2:.2f}")
            
            # Determine wedge type based on sign of t values
            if t1 > 0 and t2 > 0:
                # Both positive - up-up wedge
                wedge['type'] = 'intersecting'
                wedge['intersection_point'] = intersection
                wedge['cylinder_axis_location'] = intersection
                wedge['intersection_type'] = 'up-up'
                intersection_points.append(intersection)
                arcpy.AddMessage(f"  -> UP-UP wedge (both t positive)")
            elif t1 < 0 and t2 < 0:
                # Both negative - down-down wedge
                wedge['type'] = 'intersecting'
                wedge['intersection_point'] = intersection
                wedge['cylinder_axis_location'] = intersection
                wedge['intersection_type'] = 'down-down'
                intersection_points.append(intersection)
                arcpy.AddMessage(f"  -> DOWN-DOWN wedge (both t negative)")
            else:
                # Opposite signs - invalid wedge
                arcpy.AddMessage(f"  -> INVALID: t values have opposite signs (t1={t1:.2f}, t2={t2:.2f})")
                arcpy.AddWarning(f"No valid wedge geometry for measurements {label1} and {label2}")
        
        # Calculate plot dimensions (before including intersection points)
        x_range = profile_x_coords.max() - profile_x_coords.min()
        y_range = profile_y_coords.max() - profile_y_coords.min()
        
        # Check if all adjacent pairs have valid geometry now that intersections are calculated
        arcpy.AddMessage("\nValidating wedge coverage...")
        all_valid = True
        for wedge in wedge_data:
            if wedge['type'] not in ['parallel', 'intersecting']:
                all_valid = False
                break
        
        if not all_valid:
            arcpy.AddWarning("Stratigraphic heights cannot be calculated in the study area because some areas are not covered by valid wedges or rectangles.")
            arcpy.AddWarning("All adjacent measurement pairs must have either valid intersections or parallel geometry with consistent stratigraphic directions.")
        else:
            arcpy.AddMessage("  All adjacent pairs have valid geometry ✓")
            
            # Calculate stratigraphic heights if field was provided
            if 'strat_height' in strike_dip_data and any(h is not None for h in strike_dip_data['strat_height']):
                arcpy.AddMessage("\nCalculating stratigraphic heights...")
                self.calculate_stratigraphic_heights_in_profile(
                    strike_dip_data, wedge_data, profile_x_coords, profile_y_coords, sorted_indices)
            else:
                # No stratigraphic height field provided - assume leftmost point has height 0
                arcpy.AddMessage("\nNo stratigraphic height field provided - assuming leftmost point has stratigraphic height = 0")
                
                # Create a synthetic strat_height array with leftmost point = 0, others = None
                n_points = len(strike_dip_data['x'])
                strat_height = [None] * n_points
                
                # Find leftmost point in sorted order (first in sorted_indices)
                leftmost_idx = sorted_indices[0]
                strat_height[leftmost_idx] = 0.0
                
                strike_dip_data['strat_height'] = strat_height
                
                arcpy.AddMessage("\nCalculating stratigraphic heights...")
                self.calculate_stratigraphic_heights_in_profile(
                    strike_dip_data, wedge_data, profile_x_coords, profile_y_coords, sorted_indices)
            
            # Create output feature class with calculated heights if requested
            if output_strike_dip_fc:
                arcpy.AddMessage(f"\nCreating output strike/dip feature class with calculated heights...")
                arcpy.AddMessage(f"  Output path provided: {output_strike_dip_fc}")
                self.create_output_strike_dip_fc(strike_dip_fc, strike_dip_data, output_strike_dip_fc, label_field)
                arcpy.AddMessage(f"  Output feature class saved to: {output_strike_dip_fc}")
            else:
                arcpy.AddMessage("\nNo output feature class path provided - skipping feature class creation")
            
            # Assign DEM cells to wedges now that geometry is complete
            if dem_profile_x is not None:
                dem_wedge_assignments = self.assign_dem_cells_to_wedges(
                    dem_profile_x, dem_profile_y, strike_dip_data, wedge_data, sorted_indices)
                
                # Calculate stratigraphic heights for DEM cells
                dem_strat_heights = self.calculate_dem_stratigraphic_heights(
                    dem_profile_x, dem_profile_y, dem_wedge_assignments, strike_dip_data, wedge_data)
                
                # Create wedge assignment raster if output path provided
                if wedge_raster_output:
                    arcpy.AddMessage("\nCreating wedge assignment raster...")
                    self.create_wedge_assignment_raster(
                        dem_cell_data, dem_wedge_assignments, dem_raster, study_area, wedge_raster_output)
                    arcpy.AddMessage(f"  Wedge assignment raster saved to: {wedge_raster_output}")
                
                # Create stratigraphic height raster if output path provided
                if strat_height_raster_output:
                    arcpy.AddMessage("\nCreating stratigraphic height raster...")
                    self.create_strat_height_raster(
                        dem_cell_data, dem_strat_heights, dem_raster, study_area, strat_height_raster_output)
                    arcpy.AddMessage(f"  Stratigraphic height raster saved to: {strat_height_raster_output}")
        
        # Filter out intersection points that are too far away (>5x plot width)
        max_distance = 5 * x_range
        filtered_intersections = []
        for wedge in wedge_data:
            if wedge['intersection_point'] is not None:
                ix, iy = wedge['intersection_point']
                # Check distance from center of data
                center_x = (profile_x_coords.min() + profile_x_coords.max()) / 2
                center_y = (profile_y_coords.min() + profile_y_coords.max()) / 2
                dist = np.sqrt((ix - center_x)**2 + (iy - center_y)**2)
                
                if dist <= max_distance:
                    filtered_intersections.append((ix, iy))
                    wedge['plot_intersection'] = True
                else:
                    wedge['plot_intersection'] = False
        
        # Create the plot
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Plot DEM cells if provided AND plotting is requested, colored by wedge assignment
        if dem_profile_x is not None and dem_profile_y is not None and plot_dem_cells:
            n_dem_cells = len(dem_profile_x)
            arcpy.AddMessage(f"  Plotting DEM cells...")
            
            # Check if wedge assignments were created (may be None if validation failed)
            if dem_wedge_assignments is None:
                arcpy.AddMessage(f"  Warning: Cannot color DEM cells by wedge (wedge validation failed)")
                # Just plot all cells in gray without wedge colors
                if n_dem_cells > 10000:
                    arcpy.AddMessage(f"  Subsampling 10000 of {n_dem_cells} DEM cells for plotting...")
                    plot_indices = np.random.choice(n_dem_cells, 10000, replace=False)
                    dem_plot_x = dem_profile_x[plot_indices]
                    dem_plot_y = dem_profile_y[plot_indices]
                else:
                    arcpy.AddMessage(f"  Plotting all {n_dem_cells} DEM cells...")
                    dem_plot_x = dem_profile_x
                    dem_plot_y = dem_profile_y
                ax.scatter(dem_plot_x, dem_plot_y, c='lightgray', s=0.5, alpha=0.3, zorder=1, label='DEM Cells')
            
            else:
                # If too many cells, randomly subsample for plotting
                max_plot_cells = 10000
                if n_dem_cells > max_plot_cells:
                    arcpy.AddMessage(f"  Subsampling {max_plot_cells} of {n_dem_cells} DEM cells for plotting...")
                    plot_indices = np.random.choice(n_dem_cells, max_plot_cells, replace=False)
                    dem_plot_x = dem_profile_x[plot_indices]
                    dem_plot_y = dem_profile_y[plot_indices]
                    dem_plot_assignments = dem_wedge_assignments[plot_indices]
                else:
                    arcpy.AddMessage(f"  Plotting all {n_dem_cells} DEM cells...")
                    dem_plot_x = dem_profile_x
                    dem_plot_y = dem_profile_y
                    dem_plot_assignments = dem_wedge_assignments
                
                # Create color map for wedges
                n_wedges = len(wedge_data)
                cmap = plt.cm.get_cmap('tab20', n_wedges)
                
                # Plot points colored by wedge
                for wedge_idx in range(n_wedges):
                    mask = (dem_plot_assignments == wedge_idx)
                    if np.sum(mask) > 0:
                        ax.scatter(dem_plot_x[mask], dem_plot_y[mask], 
                                  c=[cmap(wedge_idx)], s=0.5, alpha=0.5, zorder=1,
                                  label=f'Wedge {wedge_idx}')
                
                # Plot unassigned points in gray
                mask_unassigned = (dem_plot_assignments == -1)
                if np.sum(mask_unassigned) > 0:
                    ax.scatter(dem_plot_x[mask_unassigned], dem_plot_y[mask_unassigned],
                              c='lightgray', s=0.5, alpha=0.3, zorder=1,
                              label='Unassigned')
                
                # Plot ambiguous points (multiple wedges) in black
                mask_ambiguous = (dem_plot_assignments == -2)
                if np.sum(mask_ambiguous) > 0:
                    ax.scatter(dem_plot_x[mask_ambiguous], dem_plot_y[mask_ambiguous],
                              c='black', s=1.0, alpha=0.7, zorder=1,
                              label='Ambiguous (multiple wedges)')
        
        # Calculate plot dimensions for scaling
        x_range = profile_x_coords.max() - profile_x_coords.min()
        y_range = profile_y_coords.max() - profile_y_coords.min()
        
        # Attitude line length = 3% of plot width
        line_length = 0.03 * x_range
        
        # Plot the projected points (smaller size so lines are visible)
        ax.scatter(profile_x_coords, profile_y_coords, c='red', s=25, 
                  edgecolors='black', linewidths=0.5, zorder=3, label='Measurements')
        
        # Plot projected attitude lines (blue) and stratigraphic height lines (yellow)
        for i in range(len(profile_x_coords)):
            attitude = strike_dip_data['projected_attitudes'][i]
            strat_vec = strike_dip_data['strat_height_vectors'][i]
            
            # Calculate bedding trace line endpoints centered on the point
            dx = attitude['x'] * line_length / 2
            dy = attitude['y'] * line_length / 2
            
            x1 = profile_x_coords[i] - dx
            x2 = profile_x_coords[i] + dx
            y1 = profile_y_coords[i] - dy
            y2 = profile_y_coords[i] + dy
            
            # Plot bedding trace (blue)
            ax.plot([x1, x2], [y1, y2], 'b-', linewidth=1.5, zorder=2, 
                   label='Bedding Trace' if i == 0 else '')
            
            # Calculate strat height line endpoints using stored direction
            sx1 = profile_x_coords[i]
            sy1 = profile_y_coords[i]
            sx2 = profile_x_coords[i] + strat_vec['x'] * line_length / 2
            sy2 = profile_y_coords[i] + strat_vec['y'] * line_length / 2
            
            # Plot stratigraphic height vector (yellow/gold)
            ax.plot([sx1, sx2], [sy1, sy2], color='gold', linewidth=1.5, zorder=2,
                   label='Strat Height Vector' if i == 0 else '')
        
        # Add labels for each point using the label field
        for i in range(len(profile_x_coords)):
            ax.annotate(strike_dip_data['labels'][i], (profile_x_coords[i], profile_y_coords[i]),
                       xytext=(3, 3), textcoords='offset points', fontsize=8)
        
        # Plot wedge geometry
        for wedge in wedge_data:
            idx1 = wedge['point1_idx']
            idx2 = wedge['point2_idx']
            
            # Check if these are adjacent in sorted order
            sorted_pos1 = np.where(sorted_indices == idx1)[0][0]
            sorted_pos2 = np.where(sorted_indices == idx2)[0][0]
            
            if abs(sorted_pos1 - sorted_pos2) != 1:
                continue
            
            p1_x = profile_x_coords[idx1]
            p1_y = profile_y_coords[idx1]
            p2_x = profile_x_coords[idx2]
            p2_y = profile_y_coords[idx2]
            
            if wedge['type'] == 'parallel':
                # Draw rectangle with corners at the two attitude lines
                att1 = strike_dip_data['projected_attitudes'][idx1]
                att2 = strike_dip_data['projected_attitudes'][idx2]
                
                # Rectangle corners (using attitude direction)
                dx = att1['x'] * line_length / 2
                dy = att1['y'] * line_length / 2
                
                # Four corners of rectangle
                corners_x = [p1_x - dx, p1_x + dx, p2_x + dx, p2_x - dx, p1_x - dx]
                corners_y = [p1_y - dy, p1_y + dy, p2_y + dy, p2_y - dy, p1_y - dy]
                
                ax.plot(corners_x, corners_y, 'g-', linewidth=1, alpha=0.7, zorder=1)
                
            elif wedge['type'] == 'intersecting' and wedge.get('plot_intersection', False):
                # Draw lines from points to intersection
                ix, iy = wedge['intersection_point']
                
                ax.plot([p1_x, ix], [p1_y, iy], 'k--', linewidth=0.8, alpha=0.6, zorder=1)
                ax.plot([p2_x, ix], [p2_y, iy], 'k--', linewidth=0.8, alpha=0.6, zorder=1)
                
                # Plot intersection point
                ax.plot(ix, iy, 'ko', markersize=4, zorder=2)
        
        # Format plot with padding to ensure all points are visible
        padding = 0.1  # 10% padding on each side
        
        ax.set_xlim(profile_x_coords.min() - padding * x_range, 
                    profile_x_coords.max() + padding * x_range)
        ax.set_ylim(profile_y_coords.min() - padding * y_range, 
                    profile_y_coords.max() + padding * y_range)
        
        ax.set_xlabel('Profile Distance (m)', fontsize=12)
        ax.set_ylabel('Vertical Distance (m)', fontsize=12)
        ax.set_title('Downplunge Projection of Strike/Dip Measurements', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.set_aspect('equal', adjustable='datalim')
        ax.legend(loc='best')
        
        # Invert x-axis to mirror the plot (conventional down-plunge view)
        ax.invert_xaxis()
        
        # Add scale bar
        self.add_scale_bar(ax, x_range)
        
        # Add fold axis information as text
        text_str = f'Fold Axis: {fold_trend:.1f}° / {fold_plunge:.1f}°'
        ax.text(0.02, 0.98, text_str, transform=ax.transAxes,
               fontsize=11, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        # Save to PDF
        with PdfPages(output_pdf) as pdf:
            pdf.savefig(fig, bbox_inches='tight')
            plt.close()
    
    def find_line_intersection(self, p1x, p1y, d1x, d1y, p2x, p2y, d2x, d2y):
        """
        Find intersection of two lines defined by point and direction
        Line 1: point (p1x, p1y) with direction (d1x, d1y)
        Line 2: point (p2x, p2y) with direction (d2x, d2y)
        Returns (x, y) of intersection or None if parallel
        """
        
        # Lines are: P1 + t1*D1 = P2 + t2*D2
        # Solve for t1 and t2
        # p1x + t1*d1x = p2x + t2*d2x
        # p1y + t1*d1y = p2y + t2*d2y
        
        # Rearrange: t1*d1x - t2*d2x = p2x - p1x
        #            t1*d1y - t2*d2y = p2y - p1y
        
        # Matrix form: [d1x  -d2x] [t1]   [p2x - p1x]
        #              [d1y  -d2y] [t2] = [p2y - p1y]
        
        det = d1x * (-d2y) - (-d2x) * d1y
        
        if abs(det) < 1e-10:
            # Lines are parallel
            return None
        
        dx = p2x - p1x
        dy = p2y - p1y
        
        t1 = (dx * (-d2y) - (-d2x) * dy) / det
        
        # Calculate intersection point
        ix = p1x + t1 * d1x
        iy = p1y + t1 * d1y
        
        return (ix, iy)
    
    def find_line_intersection_with_params(self, p1x, p1y, d1x, d1y, p2x, p2y, d2x, d2y):
        """
        Find intersection of two lines defined by point and direction, returning parameters
        Line 1: point (p1x, p1y) with direction (d1x, d1y)
        Line 2: point (p2x, p2y) with direction (d2x, d2y)
        Returns ((x, y), t1, t2) where t1 is parameter along line 1, t2 along line 2
        Returns (None, None, None) if parallel
        """
        
        det = d1x * (-d2y) - (-d2x) * d1y
        
        if abs(det) < 1e-10:
            # Lines are parallel
            return None, None, None
        
        dx = p2x - p1x
        dy = p2y - p1y
        
        t1 = (dx * (-d2y) - (-d2x) * dy) / det
        t2 = (d1x * dy - dx * d1y) / det
        
        # Calculate intersection point
        ix = p1x + t1 * d1x
        iy = p1y + t1 * d1y
        
        return (ix, iy), t1, t2
    
    def add_scale_bar(self, ax, x_range):
        """Add a scale bar to the plot in the lower right corner"""
        
        # Determine appropriate scale bar length based on plot range
        # Use nice round numbers: 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, etc.
        nice_numbers = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000]
        
        # Target scale bar to be about 20-30% of plot width
        target_length = x_range * 0.25
        
        # Find the closest nice number
        scale_length = min(nice_numbers, key=lambda x: abs(x - target_length))
        
        # Position scale bar in lower right corner
        # Note: x-axis is inverted, so we need to adjust accordingly
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        y_range = ylim[1] - ylim[0]
        
        # Since x-axis is inverted, xlim[0] is actually the right side
        # Scale bar position: left edge of plot (xlim[1]) plus 5% padding, then add bar length
        # Bottom edge plus 8% padding
        bar_x_start = xlim[1] + 0.05 * x_range
        bar_x_end = bar_x_start + scale_length
        bar_y = ylim[0] + 0.08 * y_range
        
        # Draw scale bar
        ax.plot([bar_x_start, bar_x_end], [bar_y, bar_y], 'k-', linewidth=2, zorder=10)
        
        # Add ticks at ends
        tick_height = 0.01 * y_range
        ax.plot([bar_x_start, bar_x_start], [bar_y - tick_height, bar_y + tick_height], 
                'k-', linewidth=2, zorder=10)
        ax.plot([bar_x_end, bar_x_end], [bar_y - tick_height, bar_y + tick_height], 
                'k-', linewidth=2, zorder=10)
        
        # Add label
        bar_x_mid = (bar_x_start + bar_x_end) / 2
        ax.text(bar_x_mid, bar_y + 3 * tick_height, f'{scale_length} m', 
                ha='center', va='bottom', fontsize=10, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

# -*- coding: utf-8 -*-
# ThreePointProblemToolbox.pyt
# ArcGIS Pro Python Toolbox for computing strike/dip from point triplets on a DEM

import arcpy
import math

arcpy.env.overwriteOutput = True


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def _ensure_fields(fc, field_specs):
    """Add fields to fc if they don't already exist.
    field_specs: list of (name, type, length_or_None)
    """
    existing = {f.name.upper() for f in arcpy.ListFields(fc)}
    for name, ftype, flen in field_specs:
        if name.upper() not in existing:
            if ftype.upper() == "TEXT" and flen:
                arcpy.management.AddField(fc, name, ftype, field_length=flen)
            else:
                arcpy.management.AddField(fc, name, ftype)


def _sample_elevation(raster, x, y):
    """Return the raster cell value at (x, y), or None if outside extent."""
    result = arcpy.GetCellValue_management(raster, f"{x} {y}")
    val = result.getOutput(0)
    if val.lower() in ("", "nodata"):
        return None
    return float(val)


def _strike_dip_from_three_points(p1, p2, p3):
    """Compute strike and dip from three (x, y, z) points.

    Returns (strike_deg, dip_deg) using the right-hand rule convention,
    or (nan, nan) if the points are collinear / coincident.
    """
    v1 = (p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2])
    v2 = (p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2])

    # Normal vector via cross product
    nx = v1[1] * v2[2] - v1[2] * v2[1]
    ny = v1[2] * v2[0] - v1[0] * v2[2]
    nz = v1[0] * v2[1] - v1[1] * v2[0]

    norm = math.sqrt(nx * nx + ny * ny + nz * nz)
    if norm == 0:
        return (float("nan"), float("nan"))

    # Force normal to point upward
    if nz < 0:
        nx, ny, nz = -nx, -ny, -nz

    # Dip: angle between the plane and horizontal
    h = math.sqrt(nx * nx + ny * ny)
    dip_deg = math.degrees(math.atan2(h, nz))

    # Down-dip azimuth, then strike by RHR (strike = dip_direction - 90)
    dip_dir = math.degrees(math.atan2(nx, ny))
    if dip_dir < 0:
        dip_dir += 360.0
    strike = (dip_dir - 90.0) % 360.0

    return (strike, dip_deg)


# ---------------------------------------------------------------------------
# Toolbox
# ---------------------------------------------------------------------------

class Toolbox(object):
    def __init__(self):
        self.label = "Three Point Problem Toolbox"
        self.alias = "threepoint"
        self.tools = [ThreePointProblem]


# ---------------------------------------------------------------------------
# Tool
# ---------------------------------------------------------------------------

class ThreePointProblem(object):
    def __init__(self):
        self.label = "Compute Strike/Dip from Point Triplets"
        self.description = (
            "Given a feature class of point triplets (ordered by OID) and an elevation "
            "raster, computes strike (0-360) and dip (0-90) for each triplet and appends "
            "a point at the centroid of each new triplet to the output feature class. "
            "Already-processed triplets are skipped by comparing the current output row "
            "count to the input point count."
        )
        self.canRunInBackground = True

    def getParameterInfo(self):
        in_pts = arcpy.Parameter(
            displayName="Point Triplets Layer",
            name="in_points",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input"
        )

        out_fc = arcpy.Parameter(
            displayName="Output Strike/Dip Layer",
            name="out_fc",
            datatype="GPFeatureLayer",
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

        return [in_pts, out_fc, in_raster]

    def isLicensed(self):
        return True

    def updateParameters(self, params):
        return

    def updateMessages(self, params):
        return

    def execute(self, params, messages):
        in_points_param = params[0].valueAsText
        out_fc_param    = params[1].valueAsText
        in_raster       = params[2].valueAsText

        # ------------------------------------------------------------------
        # Resolve layer references to catalog paths
        # ------------------------------------------------------------------
        try:
            in_points = arcpy.Describe(in_points_param).catalogPath
        except Exception:
            in_points = in_points_param

        try:
            out_fc = arcpy.Describe(out_fc_param).catalogPath
        except Exception:
            out_fc = out_fc_param

        # ------------------------------------------------------------------
        # Validate input geometry type
        # ------------------------------------------------------------------
        in_desc = arcpy.Describe(in_points)
        if getattr(in_desc, "shapeType", "").lower() != "point":
            raise arcpy.ExecuteError("The Point Triplets Layer must be a point feature class.")

        out_desc = arcpy.Describe(out_fc)
        if getattr(out_desc, "shapeType", "").lower() != "point":
            raise arcpy.ExecuteError("The Output Strike/Dip Layer must be a point feature class.")

        # ------------------------------------------------------------------
        # Ensure required output fields exist
        # ------------------------------------------------------------------
        _ensure_fields(out_fc, [
            ("STRIKE",   "DOUBLE", None),
            ("DIP",      "DOUBLE", None),
            ("SRC_OID1", "LONG",   None),
            ("SRC_OID2", "LONG",   None),
            ("SRC_OID3", "LONG",   None),
        ])

        # ------------------------------------------------------------------
        # Determine how many triplets have already been processed
        # ------------------------------------------------------------------
        n_existing = int(arcpy.management.GetCount(out_fc).getOutput(0))
        skip_points = n_existing * 3
        arcpy.AddMessage(f"Output already contains {n_existing} result(s); skipping first {skip_points} input point(s).")

        # ------------------------------------------------------------------
        # Read input points into an ordered list, skipping already-used ones
        # ------------------------------------------------------------------
        oid_field = in_desc.OIDFieldName

        # Read all points sorted by OID, skip the first skip_points,
        # then process the rest in batches of 3.
        all_points = []
        with arcpy.da.SearchCursor(
            in_points, [oid_field, "SHAPE@"],
            sql_clause=(None, f"ORDER BY {oid_field}")
        ) as cur:
            for oid, geom in cur:
                all_points.append((oid, geom))

        new_points = all_points[skip_points:]

        # Warn about incomplete trailing triplets
        remainder = len(new_points) % 3
        if remainder != 0:
            arcpy.AddWarning(
                f"{remainder} leftover point(s) at the end of the input do not form a "
                f"complete triplet and will be skipped. Add {3 - remainder} more point(s) "
                f"and re-run to process them."
            )

        # Build groups of 3
        new_groups = [
            new_points[i:i + 3]
            for i in range(0, len(new_points) - remainder, 3)
        ]

        arcpy.AddMessage(f"{len(new_groups)} new triplet(s) to process.")

        if not new_groups:
            arcpy.AddMessage("Nothing to do. No new triplets found.")
            return

        # ------------------------------------------------------------------
        # Sample elevations and compute strike/dip for each new triplet
        # ------------------------------------------------------------------
        inserted = 0

        with arcpy.da.InsertCursor(
            out_fc, ["SHAPE@XY", "STRIKE", "DIP", "SRC_OID1", "SRC_OID2", "SRC_OID3"]
        ) as icur:
            for trip in new_groups:
                oids  = [t[0] for t in trip]
                geoms = [t[1] for t in trip]

                # Sample elevation at each point
                pts_3d = []
                failed = False
                for i, (oid, geom) in enumerate(zip(oids, geoms)):
                    cx = geom.centroid.X
                    cy = geom.centroid.Y
                    z  = _sample_elevation(in_raster, cx, cy)
                    if z is None:
                        arcpy.AddWarning(
                            f"Could not sample elevation for OID {oid} at ({cx:.4f}, {cy:.4f}); "
                            f"skipping this triplet (OIDs {oids})."
                        )
                        failed = True
                        break
                    pts_3d.append((cx, cy, z))

                if failed:
                    continue

                # Compute strike and dip
                strike, dip = _strike_dip_from_three_points(*pts_3d)
                if math.isnan(strike) or math.isnan(dip):
                    arcpy.AddWarning(
                        f"Triplet OIDs {oids} are collinear or coincident; skipping."
                    )
                    continue

                # Output point at the location of the second input point
                cx_out = pts_3d[1][0]
                cy_out = pts_3d[1][1]

                icur.insertRow([(cx_out, cy_out), strike, dip, oids[0], oids[1], oids[2]])
                inserted += 1
                arcpy.AddMessage(
                    f"Triplet OIDs {oids} -> Strike={strike:.1f}, Dip={dip:.1f}"
                )

        arcpy.AddMessage(f"Done. Appended {inserted} new strike/dip point(s) to the output.")

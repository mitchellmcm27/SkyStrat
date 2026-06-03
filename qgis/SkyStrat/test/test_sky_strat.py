"""
SkyStrat plugin tests — unit and integration tests in one file.
Run with: pytest qgis/sky_strat/test/test_sky_strat.py
"""

# ---------------------------------------------------------------------------
# Stub qgis/PyQt so algorithm modules import without a running QGIS instance.
# Must happen before any sky_strat imports.
# ---------------------------------------------------------------------------
import itertools
import math
import os
import sys
import tempfile
from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
from osgeo import gdal, ogr, osr
from qgis.testing import unittest

from SkyStrat.algorithm_busk_down_plunge import BuskDownPlungeAlgorithm
from SkyStrat.algorithm_compute_strike_dip import _strike_dip_from_three_points
from SkyStrat.algorithm_project_plane import _autodetect_field, _write_circle_shapefile
from SkyStrat.shared import _sample_raster_at_xy

gdal.UseExceptions()

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


def _make_mem_raster(data, gt, nodata=None):
    """Return an in-memory GDAL dataset from a 2-D float32 array."""
    nrows, ncols = data.shape
    drv = gdal.GetDriverByName("MEM")
    ds = drv.Create("", ncols, nrows, 1, gdal.GDT_Float32)
    ds.SetGeoTransform(gt)
    band = ds.GetRasterBand(1)
    band.WriteArray(data.astype(np.float32))
    if nodata is not None:
        band.SetNoDataValue(float(nodata))
    return ds


def _make_strike_dip_data(strikes, dips, overturned=None, xs=None, ys=None, zs=None):
    """Build a minimal strike_dip_data dict for Busk method tests."""
    n = len(strikes)
    if overturned is None:
        overturned = [0] * n
    # Default: space points N-S (y direction) so they project to distinct
    # profile_x values when fold_trend=90° (EW axis, profile_x = North).
    if xs is None:
        xs = [0.0] * n
    if ys is None:
        ys = [(i - (n - 1) / 2.0) * 500.0 for i in range(n)]
    if zs is None:
        zs = [0.0] * n
    return {
        "strike": np.array(strikes, dtype=object),
        "dip": np.array(dips, dtype=object),
        "overturned": np.array(overturned, dtype=object),
        "labels": np.array([str(i) for i in range(n)], dtype=object),
        "x": np.array(xs, dtype=object),
        "y": np.array(ys, dtype=object),
        "z": np.array(zs, dtype=object),
        "strat_height": np.array([None] * n, dtype=object),
    }


# Minimal WGS-84 / UTM zone 33N WKT for shapefile tests.
_UTM33N_WKT = (
    'PROJCS["WGS 84 / UTM zone 33N",'
    'GEOGCS["WGS 84",DATUM["WGS_1984",'
    'SPHEROID["WGS 84",6378137,298.257223563]],'
    'PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],'
    'PROJECTION["Transverse_Mercator"],'
    'PARAMETER["latitude_of_origin",0],'
    'PARAMETER["central_meridian",15],'
    'PARAMETER["scale_factor",0.9996],'
    'PARAMETER["false_easting",500000],'
    'PARAMETER["false_northing",0],'
    'UNIT["metre",1]]'
)


# ===========================================================================
# 1. _sample_raster_at_xy
# ===========================================================================


class TestSampleRasterAtXY(unittest.TestCase):
    """
    3×3 raster, values 1–9.
    GeoTransform: origin (0, 3), pixel 1 m wide, 1 m tall (north-up → gt[5]=-1).
         col 0   col 1   col 2
    row 0:  1       2       3      centroid y = 2.5
    row 1:  4       5       6      centroid y = 1.5
    row 2:  7       8       9      centroid y = 0.5
    """

    DATA = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype=np.float32)
    GT = (0.0, 1.0, 0.0, 3.0, 0.0, -1.0)

    def _ds(self, data=None, nodata=None):
        return _make_mem_raster(self.DATA if data is None else data, self.GT, nodata)

    def test_centre_of_middle_cell(self):
        self.assertEqual(_sample_raster_at_xy(self._ds(), 1.5, 1.5), 5.0)

    def test_top_left_cell(self):
        self.assertEqual(_sample_raster_at_xy(self._ds(), 0.5, 2.5), 1.0)

    def test_bottom_right_cell(self):
        self.assertEqual(_sample_raster_at_xy(self._ds(), 2.5, 0.5), 9.0)

    def test_outside_x_returns_none(self):
        self.assertNone(_sample_raster_at_xy(self._ds(), 5.0, 1.5))

    def test_outside_y_returns_none(self):
        self.assertNone(_sample_raster_at_xy(self._ds(), 1.5, 5.0))

    def test_negative_coordinate_returns_none(self):
        self.assertNone(_sample_raster_at_xy(self._ds(), -0.5, 1.5))

    def test_nodata_cell_returns_none(self):
        data = self.DATA.copy()
        data[1, 1] = -9999.0
        self.assertNone(
            _sample_raster_at_xy(self._ds(data=data, nodata=-9999.0), 1.5, 1.5)
        )

    def test_non_nodata_cell_still_readable(self):
        data = self.DATA.copy()
        data[1, 1] = -9999.0
        self.assertEqual(
            _sample_raster_at_xy(self._ds(data=data, nodata=-9999.0), 0.5, 2.5), 1.0
        )

    def test_no_nodata_set_returns_value(self):
        # No NoData configured → -9999 is just a normal value.
        data = self.DATA.copy()
        data[0, 0] = -9999.0
        self.assertEqual(
            _sample_raster_at_xy(self._ds(data=data), 0.5, 2.5), pytest.approx(-9999.0)
        )

    def test_non_square_pixels(self):
        # 2 m wide, 0.5 m tall pixels; 1 row × 2 cols.
        gt = (0.0, 2.0, 0.0, 1.5, 0.0, -0.5)
        data = np.array([[10.0, 20.0]], dtype=np.float32)
        ds = _make_mem_raster(data, gt)
        self.assertEqual(_sample_raster_at_xy(ds, 1.0, 1.25), 10.0)
        self.assertEqual(_sample_raster_at_xy(ds, 3.0, 1.25), 20.0)


# ===========================================================================
# 2. _strike_dip_from_three_points
# ===========================================================================


class TestStrikeDipFromThreePoints:
    # ── geometric correctness ──────────────────────────────────────────────

    def test_horizontal_plane_gives_zero_dip(self):
        _, dip = _strike_dip_from_three_points((0, 0, 0), (1, 0, 0), (0, 1, 0))
        assert dip == pytest.approx(0.0, abs=1e-9)

    def test_45_degree_north_dipping(self):
        # Plane: origin, 1 m east (along strike), 1 m north + 1 m down.
        # Dip direction = north (000°), so RHR strike = (000−90)%360 = 270°.
        strike, dip = _strike_dip_from_three_points((0, 0, 0), (1, 0, 0), (0, 1, -1))
        assert dip == pytest.approx(45.0, abs=1e-9)
        assert strike == pytest.approx(270.0, abs=1e-9)

    def test_30_degree_south_dipping(self):
        # Dip direction = south (180°), RHR strike = (180−90)%360 = 090°.
        slope = math.tan(math.radians(30))
        strike, dip = _strike_dip_from_three_points(
            (0, 0, 0), (1, 0, 0), (0, -1, -slope)
        )
        assert dip == pytest.approx(30.0, abs=1e-9)
        assert strike == pytest.approx(90.0, abs=1e-9)

    def test_vertical_plane(self):
        # Three points in the XZ plane → dip = 90°.
        _, dip = _strike_dip_from_three_points((0, 0, 0), (1, 0, 0), (0, 0, 1))
        assert dip == pytest.approx(90.0, abs=1e-9)

    def test_rhr_dip_direction_is_90_cw_from_strike(self):
        # For any valid plane, dip_direction = (strike + 90) % 360.
        strike, _ = _strike_dip_from_three_points((0, 0, 0), (1, 0, 0), (0, 1, -1))
        dip_dir = (strike + 90.0) % 360.0
        # Expected dip direction = 000° (north, toward third point)
        assert dip_dir == pytest.approx(0.0, abs=1e-9)

    def test_round_trip_arbitrary_strike_dip(self):
        """Build 3 points from a known strike/dip and recover it exactly."""
        known_strike = 135.0
        known_dip = 25.0
        dipdir_rad = math.radians((known_strike + 90.0) % 360.0)
        strike_rad = math.radians(known_strike)
        slope = math.tan(math.radians(known_dip))
        # Along-strike unit vector (horizontal)
        p2 = (math.sin(strike_rad), math.cos(strike_rad), 0.0)
        # Along-dip unit vector (sloping downward)
        ux, uy = math.sin(dipdir_rad), math.cos(dipdir_rad)
        p3 = (ux, uy, -slope)
        s, d = _strike_dip_from_three_points((0, 0, 0), p2, p3)
        assert s == pytest.approx(known_strike, abs=1e-6)
        assert d == pytest.approx(known_dip, abs=1e-6)

    def test_reversed_point_order_gives_same_result(self):
        """Reversing p1↔p3 flips the raw normal; enforcement must correct it."""
        pts = ((0, 0, 0), (1, 0, 0), (0, 1, -1))
        s1, d1 = _strike_dip_from_three_points(*pts)
        s2, d2 = _strike_dip_from_three_points(*reversed(pts))
        assert s1 == pytest.approx(s2, abs=1e-9)
        assert d1 == pytest.approx(d2, abs=1e-9)

    def test_all_six_permutations_equal(self):
        pts = [(0, 0, 0), (10, 0, 0), (0, 10, -5)]
        results = [
            _strike_dip_from_three_points(*p) for p in itertools.permutations(pts)
        ]
        strikes = [r[0] for r in results]
        dips = [r[1] for r in results]
        assert max(strikes) - min(strikes) == pytest.approx(0.0, abs=1e-9)
        assert max(dips) - min(dips) == pytest.approx(0.0, abs=1e-9)

    def test_utm_scale_coordinates_stable(self):
        """Large UTM coordinates must not cause floating-point collapse."""
        ox, oy = 500_000.0, 5_000_000.0
        slope = math.tan(math.radians(20))
        _, dip = _strike_dip_from_three_points(
            (ox, oy, 100.0),
            (ox + 100, oy, 100.0),
            (ox, oy + 100, 100.0 - 100 * slope),
        )
        assert dip == pytest.approx(20.0, abs=1e-6)

    def test_millimetre_elevation_difference_stable(self):
        """Sub-millimetre Z differences must not produce NaN."""
        s, d = _strike_dip_from_three_points((0, 0, 0.0), (1, 0, 0.0), (0, 1, -0.001))
        assert not math.isnan(d)
        assert d > 0.0

    def test_strike_always_in_0_to_360(self):
        """Strike output must be in [0, 360) for all dip directions."""
        for dipdir_deg in range(0, 360, 15):
            r = math.radians(dipdir_deg)
            ux, uy = math.sin(r), math.cos(r)
            # Along-strike direction is perpendicular to dip direction in XY plane
            sx, sy = math.cos(r), -math.sin(r)
            slope = math.tan(math.radians(10))
            s, _ = _strike_dip_from_three_points(
                (0, 0, 0), (sx, sy, 0), (ux, uy, -slope)
            )
            assert 0.0 <= s < 360.0, (
                f"Strike {s} out of [0,360) for dipdir {dipdir_deg}"
            )

    # ── degenerate cases ───────────────────────────────────────────────────

    def test_collinear_points_return_nan(self):
        s, d = _strike_dip_from_three_points((0, 0, 0), (1, 1, 0), (2, 2, 0))
        assert math.isnan(s) and math.isnan(d)

    def test_two_coincident_points_return_nan(self):
        s, d = _strike_dip_from_three_points((0, 0, 0), (0, 0, 0), (1, 0, 0))
        assert math.isnan(s) and math.isnan(d)

    def test_all_coincident_points_return_nan(self):
        s, d = _strike_dip_from_three_points((5, 5, 5), (5, 5, 5), (5, 5, 5))
        assert math.isnan(s) and math.isnan(d)


# ===========================================================================
# 3. Busk — find_line_intersection_with_params
# ===========================================================================


class TestFindLineIntersectionWithParams:
    def _call(self, p1x, p1y, d1x, d1y, p2x, p2y, d2x, d2y):
        return _BUSK.find_line_intersection_with_params(
            p1x, p1y, d1x, d1y, p2x, p2y, d2x, d2y
        )

    def test_perpendicular_lines_exact_intersection(self):
        # Horizontal line through y=3; vertical line through x=2 → meet at (2,3).
        pt, t1, t2 = self._call(0, 3, 1, 0, 2, 0, 0, 1)
        assert pt == pytest.approx((2.0, 3.0))
        assert t1 == pytest.approx(2.0)
        assert t2 == pytest.approx(3.0)

    def test_oblique_lines_meet_at_origin(self):
        # Both lines pass through (0,0).
        pt, t1, t2 = self._call(-1, -1, 1, 1, 1, -1, -1, 1)
        assert pt == pytest.approx((0.0, 0.0), abs=1e-9)

    def test_returned_point_lies_on_both_lines(self):
        """Intersection satisfies p + t*d for both lines."""
        p1x, p1y, d1x, d1y = 1.0, 2.0, 3.0, 1.0
        p2x, p2y, d2x, d2y = 4.0, 0.0, -1.0, 2.0
        pt, t1, t2 = self._call(p1x, p1y, d1x, d1y, p2x, p2y, d2x, d2y)
        assert pt is not None
        assert pt[0] == pytest.approx(p1x + t1 * d1x, abs=1e-9)
        assert pt[1] == pytest.approx(p1y + t1 * d1y, abs=1e-9)
        assert pt[0] == pytest.approx(p2x + t2 * d2x, abs=1e-9)
        assert pt[1] == pytest.approx(p2y + t2 * d2y, abs=1e-9)

    def test_parallel_lines_return_none(self):
        pt, t1, t2 = self._call(0, 0, 1, 0, 0, 1, 1, 0)
        assert pt is None and t1 is None and t2 is None

    def test_antiparallel_lines_return_none(self):
        pt, t1, t2 = self._call(0, 0, 1, 0, 0, 1, -1, 0)
        assert pt is None and t1 is None and t2 is None

    def test_nearly_parallel_det_below_threshold_returns_none(self):
        # det = d1x*d2y - d2x*d1y = 1.0 * 1e-11 - 1.0 * 0 = 1e-11 < 1e-10
        pt, t1, t2 = self._call(0, 0, 1.0, 0.0, 1, 1, 1.0, 1e-11)
        assert pt is None

    def test_up_up_wedge_both_t_positive(self):
        # p1 + t1*d1 = intersection: start at (−1,0) going right (+1,0).
        # p2 + t2*d2 = intersection: start at (+1,0) going left (−1,0).
        # Both reach (0,0) at t=1 → t1>0, t2>0.
        pt, t1, t2 = self._call(-1, 0, 1, 0, 1, 0, -1, 0)
        assert pt == pytest.approx((0.0, 0.0), abs=1e-9)
        assert t1 > 0 and t2 > 0

    def test_down_down_wedge_both_t_negative(self):
        # Start beyond the intersection and move away from it.
        # p1=(1,0), d1=(1,0) → intersection at (0,0) requires t1=−1.
        # p2=(−1,0), d2=(−1,0) → intersection at (0,0) requires t2=−1 ...
        # Actually: p2=(−1,0), d2=(−1,0): pt = (−1 + t2*(−1), 0) = (0,0) → t2=−1.
        # But line 2 is the same as line 1 → parallel. Use non-collinear lines.
        # p1=(2,1), d1=(1,0); p2=(0,3), d2=(0,1).  Intersection at (2,3): t1=0, t2=0...
        # Better: force both t<0 by placing points past the intersection.
        # Line 1: p1=(3,0), d1=(1,0); intersection at (0,0) → t1=−3.
        # Line 2: p2=(0,2), d2=(0,1); intersection at (0,0) → t2=−2.
        pt, t1, t2 = self._call(3, 0, 1, 0, 0, 2, 0, 1)
        # Intersection: p1+t1*(1,0)=(3+t1, 0); p2+t2*(0,1)=(0, 2+t2).
        # 3+t1=0 → t1=−3; 2+t2=0 → t2=−2.
        assert pt == pytest.approx((0.0, 0.0), abs=1e-9)
        assert t1 < 0 and t2 < 0


# ===========================================================================
# 4. Busk — calculate_fold_axis
# ===========================================================================


class TestCalculateFoldAxis:
    def test_horizontal_ew_fold_axis(self):
        """Symmetric anticline (N-dipping left limb, S-dipping right) → EW axis, plunge≈0."""
        data = _make_strike_dip_data(
            strikes=[270, 270, 90, 90],
            dips=[30, 45, 30, 45],
        )
        trend, plunge = _BUSK.calculate_fold_axis(data)
        assert plunge == pytest.approx(0.0, abs=2.0)
        # Axis could be reported as 090° or 270°
        assert min(abs(trend - 90), abs(trend - 270)) == pytest.approx(0.0, abs=2.0)

    def test_plunging_fold_axis_recovery(self):
        """Synthetic data from a 030°/25° fold axis should be recovered within 2°."""
        fold_trend_true = 30.0
        fold_plunge_true = 25.0
        tr = math.radians(fold_trend_true)
        pr = math.radians(fold_plunge_true)
        axis = np.array(
            [math.cos(pr) * math.sin(tr), math.cos(pr) * math.cos(tr), math.sin(pr)]
        )

        # Build 4 poles by rotating a vector perpendicular to the axis.
        # Rodrigues rotation about `axis` through angles 0°, 45°, 90°, 135°.
        perp0 = np.array([math.cos(tr), -math.sin(tr), 0.0])
        perp0 /= np.linalg.norm(perp0)

        strikes, dips = [], []
        for angle_deg in [0, 45, 90, 135]:
            theta = math.radians(angle_deg)
            k = axis
            v = perp0
            pole = (
                v * math.cos(theta)
                + np.cross(k, v) * math.sin(theta)
                + k * float(np.dot(k, v)) * (1 - math.cos(theta))
            )
            pole /= np.linalg.norm(pole)
            pp = math.degrees(float(np.arcsin(np.clip(pole[2], -1, 1))))
            pt = math.degrees(math.atan2(float(pole[0]), float(pole[1])))
            dip = 90.0 - pp
            strike = (pt + 90.0) % 360.0
            strikes.append(strike)
            dips.append(max(0.0, min(90.0, dip)))

        data = _make_strike_dip_data(strikes=strikes, dips=dips)
        trend, plunge = _BUSK.calculate_fold_axis(data)

        plunge_ok = abs(plunge - fold_plunge_true) < 2.0
        trend_ok = (
            min(
                abs(trend - fold_trend_true), abs(trend - (fold_trend_true + 180) % 360)
            )
            < 2.0
        )
        assert plunge_ok and trend_ok

    def test_single_measurement_does_not_crash(self):
        data = _make_strike_dip_data(strikes=[90], dips=[30])
        trend, plunge = _BUSK.calculate_fold_axis(data)
        assert 0.0 <= plunge <= 90.0

    def test_identical_measurements_do_not_crash(self):
        data = _make_strike_dip_data(strikes=[45, 45, 45], dips=[20, 20, 20])
        trend, plunge = _BUSK.calculate_fold_axis(data)
        assert 0.0 <= plunge <= 90.0


# ===========================================================================
# 5. Busk — calculate_projected_attitudes
# ===========================================================================


class TestCalculateProjectedAttitudes:
    """All tests use fold_trend=90° (EW), fold_plunge=0° (horizontal)."""

    FOLD_TREND = 90.0
    FOLD_PLUNGE = 0.0

    def _run(self, strikes, dips, overturned=None):
        data = _make_strike_dip_data(strikes, dips, overturned=overturned)
        return _BUSK.calculate_projected_attitudes(
            data, self.FOLD_TREND, self.FOLD_PLUNGE
        )

    def test_projected_attitudes_are_unit_length(self):
        data = self._run([270, 90, 270, 90], [30, 30, 45, 45])
        for att in data["projected_attitudes"]:
            length = math.sqrt(att["x"] ** 2 + att["y"] ** 2)
            assert length == pytest.approx(1.0, abs=1e-9)

    def test_strat_vector_perpendicular_to_attitude(self):
        data = self._run([270, 90, 270, 90], [30, 30, 45, 45])
        for att, sv in zip(data["projected_attitudes"], data["strat_height_vectors"]):
            dot = att["x"] * sv["x"] + att["y"] * sv["y"]
            assert dot == pytest.approx(0.0, abs=1e-9)

    def test_strat_vector_is_unit_length(self):
        data = self._run([270, 90], [30, 45])
        for sv in data["strat_height_vectors"]:
            length = math.sqrt(sv["x"] ** 2 + sv["y"] ** 2)
            assert length == pytest.approx(1.0, abs=1e-9)

    def test_upright_beds_strat_vector_points_up(self):
        data = self._run([270, 90], [30, 30], overturned=[0, 0])
        for sv in data["strat_height_vectors"]:
            assert sv["y"] > 0.0

    def test_overturned_beds_strat_vector_points_down(self):
        data = self._run([270, 90], [30, 30], overturned=[1, 1])
        for sv in data["strat_height_vectors"]:
            assert sv["y"] < 0.0

    def test_overturned_flag_changes_projected_attitude(self):
        """Same strike/dip with overturned=1 must give a different attitude direction."""
        d_up = self._run([90], [30], overturned=[0])
        d_ov = self._run([90], [30], overturned=[1])
        att_up = d_up["projected_attitudes"][0]
        att_ov = d_ov["projected_attitudes"][0]
        dot = att_up["x"] * att_ov["x"] + att_up["y"] * att_ov["y"]
        assert dot != pytest.approx(1.0, abs=0.01)

    def test_profile_coords_have_correct_length(self):
        data = self._run([270, 90, 270], [30, 30, 45])
        assert len(data["profile_x"]) == 3
        assert len(data["profile_y"]) == 3

    def test_points_separated_in_profile_x(self):
        """N-S spaced points (default) must produce distinct profile_x values."""
        data = self._run([270, 270], [30, 30])
        assert data["profile_x"][0] != pytest.approx(data["profile_x"][1], abs=1.0)


# ===========================================================================
# 6. Busk — analyze_wedges
# ===========================================================================


class TestAnalyzeWedges:
    def _run_wedge(self, strikes, dips, overturned=None):
        data = _make_strike_dip_data(strikes, dips, overturned=overturned)
        data = _BUSK.calculate_projected_attitudes(data, 90.0, 0.0)
        sorted_indices = np.argsort(data["profile_x"])
        wedge_data = _BUSK.analyze_wedges(data, sorted_indices)
        return wedge_data, data, sorted_indices

    def test_identical_upright_attitudes_parallel(self):
        """Two identical upright measurements → 'parallel'."""
        wedge_data, _, _ = self._run_wedge([270, 270], [30, 30])
        assert len(wedge_data) == 1
        assert wedge_data[0]["type"] == "parallel"

    def test_identical_attitudes_opposite_younging_invalid_parallel(self):
        """One upright, one overturned (same geometry) → 'invalid_parallel'."""
        wedge_data, _, _ = self._run_wedge([270, 270], [30, 30], overturned=[0, 1])
        assert len(wedge_data) == 1
        assert wedge_data[0]["type"] == "invalid_parallel"

    def test_different_attitudes_type_none(self):
        """Opposite fold limbs → angle >> 0.01° → type remains None."""
        wedge_data, _, _ = self._run_wedge([270, 90], [45, 45])
        assert len(wedge_data) == 1
        assert wedge_data[0]["type"] is None

    def test_three_points_produce_two_wedges(self):
        wedge_data, _, _ = self._run_wedge([270, 90, 270], [30, 30, 45])
        assert len(wedge_data) == 2

    def test_wedge_indices_reference_sorted_adjacent_pairs(self):
        wedge_data, data, sorted_indices = self._run_wedge([270, 90, 270], [30, 30, 45])
        for i, wedge in enumerate(wedge_data):
            idx1 = wedge["point1_idx"]
            idx2 = wedge["point2_idx"]
            # idx1 and idx2 must be adjacent in the sorted order
            pos1 = int(np.where(sorted_indices == idx1)[0][0])
            pos2 = int(np.where(sorted_indices == idx2)[0][0])
            assert abs(pos1 - pos2) == 1


# ===========================================================================
# 7. Busk — assign_dem_cells_to_wedges + calculate_dem_stratigraphic_heights
# ===========================================================================


class TestWedgeAssignmentAndHeights:
    """Uses a parallel-wedge setup (two identical measurements) for isolation."""

    def _build_parallel_setup(self):
        """Return (data, wedge_data, sorted_indices) for a guaranteed parallel wedge."""
        data = _make_strike_dip_data([270, 270], [30, 30])
        data = _BUSK.calculate_projected_attitudes(data, 90.0, 0.0)
        sorted_indices = np.argsort(data["profile_x"])
        wedge_data = _BUSK.analyze_wedges(data, sorted_indices)
        assert wedge_data[0]["type"] == "parallel"
        return data, wedge_data, sorted_indices

    def test_midpoint_cell_is_assigned(self):
        data, wedge_data, sorted_indices = self._build_parallel_setup()
        px = data["profile_x"]
        py = data["profile_y"]
        mid_x = float((px[0] + px[1]) / 2.0)
        mid_y = float((py[0] + py[1]) / 2.0)

        assignments = _BUSK.assign_dem_cells_to_wedges(
            np.array([mid_x]),
            np.array([mid_y]),
            data,
            wedge_data,
            sorted_indices,
            _FB,
        )
        assert assignments[0] >= 0

    def test_far_cell_is_unassigned(self):
        data, wedge_data, sorted_indices = self._build_parallel_setup()
        assignments = _BUSK.assign_dem_cells_to_wedges(
            np.array([99999.0]),
            np.array([99999.0]),
            data,
            wedge_data,
            sorted_indices,
            _FB,
        )
        assert assignments[0] == -1

    def test_multiple_cells_some_inside_some_outside(self):
        data, wedge_data, sorted_indices = self._build_parallel_setup()
        px = data["profile_x"]
        py = data["profile_y"]
        mid_x = float((px[0] + px[1]) / 2.0)
        mid_y = float((py[0] + py[1]) / 2.0)

        dem_px = np.array([mid_x, 99999.0])
        dem_py = np.array([mid_y, 99999.0])
        assignments = _BUSK.assign_dem_cells_to_wedges(
            dem_px,
            dem_py,
            data,
            wedge_data,
            sorted_indices,
            _FB,
        )
        assert assignments[0] >= 0
        assert assignments[1] == -1

    def test_parallel_height_formula(self):
        """Height at a profile point = left_height + projection onto strat vector."""
        data, wedge_data, sorted_indices = self._build_parallel_setup()
        left_idx = int(sorted_indices[0])
        left_height = 100.0
        data["calculated_strat_height"] = [None, None]
        data["calculated_strat_height"][left_idx] = left_height

        px = data["profile_x"]
        py = data["profile_y"]
        left_x = float(px[left_idx])
        left_y = float(py[left_idx])
        right_x = float(px[int(sorted_indices[1])])
        right_y = float(py[int(sorted_indices[1])])
        mid_x = (left_x + right_x) / 2.0
        mid_y = (left_y + right_y) / 2.0

        dem_px = np.array([mid_x])
        dem_py = np.array([mid_y])

        assignments = _BUSK.assign_dem_cells_to_wedges(
            dem_px,
            dem_py,
            data,
            wedge_data,
            sorted_indices,
            _FB,
        )
        heights = _BUSK.calculate_dem_stratigraphic_heights(
            dem_px,
            dem_py,
            assignments,
            data,
            wedge_data,
            _FB,
        )

        # Manual formula (mirrors the parallel branch in the source)
        sv = data["strat_height_vectors"][left_idx]
        dm_tx = mid_x - left_x
        dm_ty = mid_y - left_y
        angle = math.atan2(sv["x"], sv["y"])
        cos_a = math.cos(-angle)
        sin_a = math.sin(-angle)
        dem_rot_y = sin_a * dm_tx + cos_a * dm_ty
        expected = dem_rot_y + left_height

        assert not np.isnan(heights[0])
        assert heights[0] == pytest.approx(expected, abs=1e-6)

    def test_unassigned_cell_height_is_nan(self):
        data, wedge_data, sorted_indices = self._build_parallel_setup()
        left_idx = int(sorted_indices[0])
        data["calculated_strat_height"] = [None, None]
        data["calculated_strat_height"][left_idx] = 100.0

        # Cell that won't be inside any wedge
        assignments = np.array([-1], dtype=int)
        heights = _BUSK.calculate_dem_stratigraphic_heights(
            np.array([99999.0]),
            np.array([99999.0]),
            assignments,
            data,
            wedge_data,
            _FB,
        )
        assert np.isnan(heights[0])

    def test_missing_calculated_heights_returns_nan_array(self):
        """If 'calculated_strat_height' is absent, all heights are NaN."""
        data, wedge_data, sorted_indices = self._build_parallel_setup()
        # Do NOT set calculated_strat_height
        px = data["profile_x"]
        py = data["profile_y"]
        mid_x = float((px[0] + px[1]) / 2.0)
        mid_y = float((py[0] + py[1]) / 2.0)
        assignments = _BUSK.assign_dem_cells_to_wedges(
            np.array([mid_x]),
            np.array([mid_y]),
            data,
            wedge_data,
            sorted_indices,
            _FB,
        )
        heights = _BUSK.calculate_dem_stratigraphic_heights(
            np.array([mid_x]),
            np.array([mid_y]),
            assignments,
            data,
            wedge_data,
            _FB,
        )
        assert np.all(np.isnan(heights))


# ===========================================================================
# 8. _autodetect_field
# ===========================================================================


class _MockField:
    def __init__(self, name):
        self._name = name

    def name(self):
        return self._name


class _MockSource:
    def __init__(self, field_names):
        self._fields = [_MockField(n) for n in field_names]

    def fields(self):
        return self._fields


class TestAutodetectField:
    def test_provided_field_returned_unchanged(self):
        src = _MockSource(["foo", "strike"])
        assert _autodetect_field(src, "strike", ["strike"]) == "strike"

    def test_candidate_found_case_insensitive(self):
        src = _MockSource(["DIP", "STRIKE"])
        assert _autodetect_field(src, None, ["dip", "dip_deg"]) == "DIP"

    def test_second_candidate_found_when_first_absent(self):
        src = _MockSource(["strike_deg"])
        assert (
            _autodetect_field(src, None, ["strike", "strike_deg", "strk"])
            == "strike_deg"
        )

    def test_no_match_returns_none(self):
        src = _MockSource(["elevation", "slope"])
        assert _autodetect_field(src, None, ["strike", "strike_deg", "strk"]) is None

    def test_provided_takes_precedence_over_candidates(self):
        src = _MockSource(["strike", "dip"])
        # 'dip' would match candidate 'dip', but provided='strike' wins.
        assert _autodetect_field(src, "strike", ["dip"]) == "strike"

    def test_empty_field_list_returns_none(self):
        src = _MockSource([])
        assert _autodetect_field(src, None, ["strike"]) is None

    def test_provided_empty_string_falls_through_to_candidates(self):
        # provided='' is falsy → treated as not provided.
        src = _MockSource(["DIP"])
        assert _autodetect_field(src, "", ["dip"]) == "DIP"


# ===========================================================================
# 9. _write_circle_shapefile
# ===========================================================================


class TestWriteCircleShapefile:
    def _make_circle(self, cx=0.0, cy=0.0, radius=100.0):
        path = tempfile.mktemp(suffix=".shp")
        geom = _write_circle_shapefile(cx, cy, radius, _UTM33N_WKT, path)
        return path, geom, radius

    def test_creates_exactly_one_polygon_feature(self):
        path, _, _ = self._make_circle()
        ds = ogr.Open(path)
        assert ds is not None
        lyr = ds.GetLayer(0)
        assert lyr.GetFeatureCount() == 1
        feat = lyr.GetNextFeature()
        assert feat.GetGeometryRef().GetGeometryType() == ogr.wkbPolygon
        ds = None

    def test_area_within_half_percent_of_pi_r_squared(self):
        path, _, radius = self._make_circle(radius=500.0)
        ds = ogr.Open(path)
        feat = ds.GetLayer(0).GetNextFeature()
        area = feat.GetGeometryRef().GetArea()
        expected = math.pi * radius**2
        assert area == pytest.approx(expected, rel=0.005)
        ds = None

    def test_envelope_centred_at_specified_point(self):
        cx, cy = 1000.0, 2000.0
        _, geom, radius = self._make_circle(cx=cx, cy=cy, radius=50.0)
        env = geom.GetEnvelope()  # (minX, maxX, minY, maxY)
        assert (env[0] + env[1]) / 2 == pytest.approx(cx, abs=1.0)
        assert (env[2] + env[3]) / 2 == pytest.approx(cy, abs=1.0)

    def test_crs_matches_input_wkt(self):
        path, _, _ = self._make_circle()
        ds = ogr.Open(path)
        srs_file = ds.GetLayer(0).GetSpatialRef()
        srs_ref = osr.SpatialReference()
        srs_ref.ImportFromWkt(_UTM33N_WKT)
        assert srs_file.GetUTMZone() == srs_ref.GetUTMZone()
        ds = None

    def test_returned_geometry_area_matches_shapefile(self):
        path, returned_geom, _ = self._make_circle()
        ds = ogr.Open(path)
        file_geom = ds.GetLayer(0).GetNextFeature().GetGeometryRef()
        assert returned_geom.GetArea() == pytest.approx(file_geom.GetArea(), rel=1e-6)
        ds = None

    def test_overwrites_existing_file(self):
        path, _, _ = self._make_circle(radius=10.0)
        # Write a second (larger) circle to the same path.
        _write_circle_shapefile(0.0, 0.0, 200.0, _UTM33N_WKT, path)
        ds = ogr.Open(path)
        feat = ds.GetLayer(0).GetNextFeature()
        area = feat.GetGeometryRef().GetArea()
        ds = None
        assert area == pytest.approx(math.pi * 200.0**2, rel=0.005)


def run_all():
    """Default function that is called by the runner if nothing else is specified"""
    suite = unittest.TestSuite()
    suite.addTests(unittest.makeSuite(TestSampleRasterAtXY, "test"))
    unittest.TextTestRunner(verbosity=3, stream=sys.stdout).run(suite)

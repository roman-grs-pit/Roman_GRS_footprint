"""Tests for sky coordinate utilities."""

import numpy as np

from rstgrs_footprint.sky_coords import generate_randoms, tangent_plane


def test_tangent_plane_center_maps_to_origin():
    """The pointing center should project to (0, 0)."""
    x, y = tangent_plane(10.0, -20.0, 10.0, -20.0, pointing_pa=0.0, focal_pa=0.0)
    assert np.isclose(x, 0.0)
    assert np.isclose(y, 0.0)


def test_tangent_plane_small_offsets_follow_expected_axes():
    """Small RA/Dec offsets should map to the expected focal-plane directions."""
    x_east, y_east = tangent_plane(11.0, -20.0, 10.0, -20.0, pointing_pa=0.0, focal_pa=0.0)
    x_north, y_north = tangent_plane(10.0, -19.0, 10.0, -20.0, pointing_pa=0.0, focal_pa=0.0)

    assert x_east > 0.0
    assert y_east < 0.0
    assert y_north > 0.0
    assert np.isclose(x_north, 0.0, atol=1e-12)


def test_tangent_plane_accepts_array_inputs():
    """The projection should work on vector inputs and preserve shape."""
    ra = np.array([10.0, 11.0, 10.0])
    dec = np.array([-20.0, -20.0, -19.0])

    x, y = tangent_plane(ra, dec, 10.0, -20.0, pointing_pa=0.0, focal_pa=0.0)

    assert x.shape == ra.shape
    assert y.shape == dec.shape
    assert np.isclose(x[0], 0.0)
    assert np.isclose(y[0], 0.0)
    assert x[1] > 0.0
    assert y[2] > 0.0


def test_tangent_plane_applies_rotation():
    """A 90 degree PA rotation should swap the projected east/north axes."""
    x_unrotated, y_unrotated = tangent_plane(11.0, -20.0, 10.0, -20.0, pointing_pa=0.0, focal_pa=0.0)
    x_rotated, y_rotated = tangent_plane(11.0, -20.0, 10.0, -20.0, pointing_pa=90.0, focal_pa=0.0)

    assert np.isclose(x_rotated, -y_unrotated, atol=1e-12)
    assert np.isclose(y_rotated, x_unrotated, atol=1e-12)


def test_generate_randoms_is_reproducible_and_within_bounds():
    """Random coordinate generation should be reproducible and bounded."""
    ra1, dec1 = generate_randoms(nran=5, ra_bounds=(10, 20), dec_bounds=(-30, 30), random_seed=7)
    ra2, dec2 = generate_randoms(nran=5, ra_bounds=(10, 20), dec_bounds=(-30, 30), random_seed=7)

    assert np.allclose(ra1, ra2)
    assert np.allclose(dec1, dec2)
    assert np.all((ra1 >= 10) & (ra1 <= 20))
    assert np.all((dec1 >= -30) & (dec1 <= 30))

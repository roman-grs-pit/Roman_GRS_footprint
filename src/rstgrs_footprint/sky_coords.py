import numpy as np

def tangent_plane(ra, dec, pointing_ra, pointing_dec, pointing_pa=0.0, focal_pa=-60.0):
    """Project (ra, dec) onto a tangent plane at a reference pointing.

    Uses the gnomonic (tangent-plane) projection with an optional
    position-angle rotation and focal plane rotation.

    Parameters
    ----------
    ra, dec : array_like
        Right ascension and declination of the source(s) in degrees.
    pointing_ra, pointing_dec : float
        Reference point (field center) in degrees.
    pointing_pa : float, optional
        Position angle in degrees, measured East of North (default 0).
    focal_pa : float, optional
        Physical orientation of the focal plane coordinate system in
        degrees, measured East of North relative to the telescope
        (default -60).  For example, -60 means the focal plane axes
        are rotated 60° West of the telescope's V3 axis.  The total
        sky orientation is PA + focal_pa; the code applies the inverse
        rotation -(PA + focal_pa) to transform sky coordinates into
        the focal plane frame.

    Returns
    -------
    x, y : ndarray
        Focal-plane coordinates in degrees.  With PA=0 and focal_pa=0,
        x points East and y points North.
    """
    # Convert everything to radians
    ra_rad = np.deg2rad(np.asarray(ra, dtype=float))
    dec_rad = np.deg2rad(np.asarray(dec, dtype=float))
    ra0 = np.deg2rad(pointing_ra)
    dec0 = np.deg2rad(pointing_dec)
    pa = np.deg2rad(-(pointing_pa + focal_pa))

    # ------------------------------------------------------------------
    # Step 1 + 2a: Shift RA so the reference point is at RA = 0,
    # then convert to 3D Cartesian coordinates on the unit sphere.
    # ------------------------------------------------------------------
    da = ra_rad - ra0  # relative RA
    cos_dec = np.cos(dec_rad)
    x = cos_dec * np.cos(da)
    y = cos_dec * np.sin(da)
    z = np.sin(dec_rad)

    # ------------------------------------------------------------------
    # Step 2b: Rotate about the y-axis by angle = dec0 - pi/2
    # to bring the reference point (which is now at dec=dec0, RA=0)
    # to the North Pole.
    #
    # The rotation angle is  theta_y = dec0 - pi/2, so:
    #   cos(theta_y) =  cos(dec0 - pi/2) =  sin(dec0)
    #   sin(theta_y) =  sin(dec0 - pi/2) = -cos(dec0)
    # ------------------------------------------------------------------
    sin_dec0 = np.sin(dec0)
    cos_dec0 = np.cos(dec0)

    # Apply R_y(dec0 - pi/2):
    #   x' =  cos(theta_y) * x + sin(theta_y) * z =  sin(dec0)*x - cos(dec0)*z
    #   y' =  y                                     (unchanged)
    #   z' = -sin(theta_y) * x + cos(theta_y) * z =  cos(dec0)*x + sin(dec0)*z
    xr = sin_dec0 * x - cos_dec0 * z
    yr = y
    zr = cos_dec0 * x + sin_dec0 * z

    # ------------------------------------------------------------------
    # Step 3: Gnomonic projection onto the tangent plane at the pole.
    # The y-axis rotation leaves yr pointing East, and xr points -Dec.
    # So to get xi=East, eta=North:
    #   xi  =  yr / zr   (East)
    #   eta = -xr / zr   (North)
    # ------------------------------------------------------------------
    xi = yr / zr
    eta = -xr / zr

    # ------------------------------------------------------------------
    # Step 4: Rotate by -(PA + focal_pa).
    # The FPA is oriented at (PA + focal_pa) on the sky. To transform
    # from sky to FPA, we apply the inverse: -(PA + focal_pa).
    # The matrix [cos, sin; -sin, cos] rotates N toward E for positive
    # angles.
    # ------------------------------------------------------------------
    cos_pa = np.cos(pa)
    sin_pa = np.sin(pa)
    x_out = cos_pa * xi + sin_pa * eta
    y_out = -sin_pa * xi + cos_pa * eta

    # Convert from radians to degrees
    return np.rad2deg(x_out), np.rad2deg(y_out)

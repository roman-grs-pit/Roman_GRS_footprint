Roman GRS footprint utilities
=============================

This package contains the coordinate and detector-footprint tools used to
translate between sky positions and the Roman focal-plane geometry.

The two most important modules are ``rstgrs_footprint.sky_coords`` and
``rstgrs_footprint.trace_on_det``. The first converts celestial coordinates into
local tangent-plane/focal-plane coordinates and generates sample sky positions.
The second checks whether a focal-plane point falls within the active detector
footprint by tracing it through the Roman optical model.

Overview
--------

``sky_coords.py``
~~~~~~~~~~~~~~~~~

The functions in ``sky_coords.py`` project RA/Dec values onto a local tangent
plane centered on a chosen pointing. This is the step that turns a sky position
into a coordinate system suitable for evaluating Roman field geometry.

``tangent_plane()`` performs a gnomonic projection from sky coordinates to
focal-plane coordinates, with optional rotation by the pointing position angle
and the focal-plane orientation. ``generate_randoms()`` creates reproducible
random RA/Dec samples within a requested sky region for testing and simulation
work.

A typical workflow is to choose a pointing center, transform sky coordinates into
local focal-plane coordinates, and then examine whether those coordinates fall
inside the instrument footprint.

.. code-block:: python

  from rstgrs_footprint import sky_coords

  x, y = sky_coords.tangent_plane(
      ra=[10.0, 11.0],
      dec=[-20.0, -20.0],
      pointing_ra=10.0,
      pointing_dec=-20.0,
      pointing_pa=0.0,
      focal_pa=-60.0,
  )

  ra_rand, dec_rand = sky_coords.generate_randoms(
      nran=100,
      ra_bounds=(0.0, 360.0),
      dec_bounds=(-90.0, 90.0),
      random_seed=42,
  )

``trace_on_det.py``
~~~~~~~~~~~~~~~~~~~

The functions in ``trace_on_det.py`` evaluate whether a focal-plane point lies
within an active detector footprint and therefore can be traced through the
instrument optics.

``test_foot()`` takes focal-plane coordinates, detector information, and
wavelength limits, then asks the Roman optical model to compute the detector
trace and returns ``1`` if the trace stays inside the detector bounds and ``0``
otherwise. This is used after converting sky coordinates into focal-plane
coordinates to answer the practical question: "Is this source on a valid
detector region for the requested wavelength range?"

.. code-block:: python

  from rstgrs_footprint import sky_coords, trace_on_det

  x, y = sky_coords.tangent_plane(ra=12.5, dec=-19.5, pointing_ra=10.0, pointing_dec=-20.0)
  on_detector = trace_on_det.test_foot(
      x,
      y,
      det=1,
      min_lam_4foot=1.0,
      max_lam_4foot=1.93,
  )

Together, these modules provide the core geometry chain needed to sample sky
positions, transform them into the Roman focal plane, and test whether they land
within the detector footprint.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

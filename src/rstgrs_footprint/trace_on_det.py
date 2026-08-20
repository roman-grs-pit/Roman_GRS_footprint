from roman_gdps_optical_model import RomanOpticalModel
import numpy as np
import os
optmod_config = os.getenv('ROMAN_GDPS_OPTICAL_MODEL_CONFIG', None)
if optmod_config is None:
    raise ValueError(
        "ROMAN_GDPS_OPTICAL_MODEL_CONFIG environment variable is not set.")
optmod = RomanOpticalModel(config_file=optmod_config)


def test_foot(xfpa, yfpa, min_pix=0, max_pix=4088, det=1, min_lam_4foot=1., max_lam_4foot=1.93):
    '''
    Tests if a point in the focal plane is within the detector bounds.

    Parameters
    ----------
    xfpa,yapa : float
        x and y position in the focal plane in degrees
    min_pix,max_pix : int
        the minimum and maximum pixel values for the detector, used to define
        the detector bounds in pixels
    det : int
        the SCA detector number
    min_lam_4foot : float
        the minimum wavelength required to be considered
    max_lam_4foot : float
        the maximum wavelength required to be considered

    Returns
    -------
    1 if the point is within the detector bounds, 0 otherwise
    '''
    optmod.wl_grid = np.array([min_lam_4foot, max_lam_4foot])
    test = optmod._get_beam_trace(xfpa, yfpa, det, width=1)
    tracex = test['trace_sca_x'][0]
    tracey = test['trace_sca_y'][0]
    if np.min(tracex) >= min_pix and np.max(tracex) < max_pix and np.min(tracey) >= min_pix and np.max(tracey) < max_pix:
        return 1
    else:
        return 0

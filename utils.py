import numpy as np
from scipy.signal import savgol_filter


def evaluate_hv_snr(x_fit, y_fit, g3, g4):
    """
    Estimate signal-to-noise ratio (SNR) for the HV component.

    Parameters:
    - x_fit : np.ndarray
        Wavelength array within the fitting window.
    - y_fit : np.ndarray
        Normalized flux values.
    - g3, g4 : np.ndarray
        Gaussian model components representing the HV feature.

    Returns:
    - snr : float
        Estimated SNR for the HV feature.
    """
    noise = np.std(y_fit - savgol_filter(y_fit, window_length=15, polyorder=2))
    signal = np.max(g3 + g4)
    return signal / noise if noise > 0 else 0



def trapezoid_area(y, x):
    """
    Compute area under a curve using trapezoidal integration.

    Parameters:
    - y : np.ndarray
        Y-values (e.g., 1 - model).
    - x : np.ndarray
        X-values (wavelengths).

    Returns:
    - area : float
        Integrated area, e.g., pseudo-equivalent width.
    """
    from numpy import trapezoid
    return trapezoid(y, x)



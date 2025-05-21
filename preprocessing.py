
import numpy as np
from scipy.ndimage import median_filter

def clean_spectrum(wave, flux, numsig=2.0, grow=1, niter=10, window=20):
    """
    Remove outliers from a spectrum using sigma-clipping and median filtering.

    Parameters:
    - wave : np.ndarray
        Wavelength or velocity array.
    - flux : np.ndarray
        Flux values to be cleaned.
    - numsig : float
        Sigma threshold for identifying outliers.
    - grow : int
        Number of neighbor pixels on each side to flag as bad.
    - niter : int
        Number of sigma-clipping iterations.
    - window : int
        Width of the smoothing kernel.

    Returns:
    - cleaned_flux : np.ndarray
        Flux array with outliers replaced by running median.
    """
    cleaned = flux.copy()
    n = len(flux)
    for _ in range(niter):
        smoothed = median_filter(cleaned, size=window, mode='nearest')
        std = np.array([np.std(cleaned[max(0, i - 10):min(n, i + 11)]) for i in range(n)])
        outliers = (cleaned > smoothed + numsig * std) | (cleaned < smoothed - numsig * std)
        if not np.any(outliers):
            break
        for idx in np.where(outliers)[0]:
            for offset in range(-grow, grow + 1):
                if 0 <= idx + offset < n:
                    cleaned[idx + offset] = smoothed[idx + offset]
    return cleaned


def fit_continuum(x, y, blue, red):
    """
    Fit linear continuum using regions on either side of the feature.

    Parameters:
    - x : np.ndarray
        Wavelength array.
    - y : np.ndarray
        Flux array.
    - blue, red : float
        Feature boundaries defining the fitting window.

    Returns:
    - continuum : np.ndarray
        Fitted linear continuum over the full x range.
    """
    cont_mask = ((x >= blue - 30) & (x <= blue)) | ((x >= red) & (x <= red + 30))
    m, b = np.polyfit(x[cont_mask], y[cont_mask], 1)
    return m * x + b

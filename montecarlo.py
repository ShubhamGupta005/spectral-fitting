
import numpy as np
from scipy.signal import savgol_filter

from preprocessing import clean_spectrum, fit_continuum
from utils import evaluate_hv_snr

def run_monte_carlo(spec, redshift, line_type="Si", n_trials=10):
    """
    Purpose:
    • Estimate uncertainties (errors) in the fitted spectral feature parameters (pEW, velocity, FWHM) using Monte Carlo resampling.
    • Adds synthetic noise to the original spectrum and refits multiple times to assess parameter variability.

    Parameter Descriptions:
    • spec : dict
        Dictionary containing 'wave' and 'flux' arrays of the observed spectrum.
    • redshift : float
        Redshift of the source for rest-frame wavelength correction.
    • line_type : str, optional
        Type of spectral line to fit ("Si" or "CaHK"). Determines rest wavelengths and fitting window. Default is "Si".
    • n_trials : int, optional
        Number of Monte Carlo trials (refits with synthetic noise). Default is 10.

    Returns:
    • result_dict : dict
        Dictionary of standard deviations for each parameter across valid trials:
        - 'pew_pv_err', 'pew_hv_err', 'pew_total_err': uncertainties in pseudo-equivalent width (Å)
        - 'vel_pv_err', 'vel_hv_err': uncertainties in velocity (km/s)
        - 'fwhm_pv_err', 'fwhm_hv_err': uncertainties in FWHM (km/s)
        - 'n_valid': number of successful trials (those not rejected due to failed fit or NaN)

    Physics Concept:
    • Propagates measurement noise into parameter uncertainties via forward-modeling.
    • Assumes the dominant source of uncertainty is Gaussian random noise in the flux.
    • This statistical approach mirrors bootstrapping or posterior sampling in error estimation, reflecting the robustness of the fit under noise.
    """
    if line_type.lower() == "si":
        from fitting import fit_siII 
        rest1, rest2 = 6347, 6371
        window, blue, red = (6000, 6500), 6050, 6275
        fit_func = fit_siII
    elif line_type.lower() == "cahk":
        from fitting import fit_cahk 
        rest1, rest2 = 3934.777, 3969.591
        window, blue, red = (3550, 4050), 3600, 4000
        fit_func = fit_cahk
    else:
        raise ValueError("line_type must be 'Si' or 'CaHK'")

    opz = 1 + redshift
    wave_rest = spec['wave'] / opz
    flux = np.array(spec['flux'])
    mask = (wave_rest >= window[0]) & (wave_rest <= window[1])
    x = wave_rest[mask]
    y_raw = flux[mask]
    resolution = np.median(np.diff(x))
    win = int(round(20.0 / resolution))
    y_clean = clean_spectrum(x, y_raw, window=win)
    y_norm = y_clean / fit_continuum(x, y_clean, blue, red)
    x_fit = x[(x > blue) & (x < red)] 
    noise = np.std(y_norm - savgol_filter(y_norm, window_length=15, polyorder=2))

    # Run MC trials
    pew_pv_list, pew_hv_list,pew_total_list = [],[],[]
    vel_pv_list, vel_hv_list = [], []
    fwhm_pv_list, fwhm_hv_list = [], []

    # Redefine x_fit to full window before blue/red masking
    x_full = wave_rest[mask]
    y_norm_full = y_clean / fit_continuum(x_full, y_clean, blue, red)

    # Add noise to full spectrum (not just narrow fit region)
    for trial in range(n_trials):
        try:
            y_noisy = y_norm_full + np.random.normal(0, noise, size=len(y_norm_full))
            spec_noisy = {"wave": x_full * (1 + redshift), "flux": y_noisy}
            result = fit_func(spec_noisy, redshift=redshift, use_snr_filter=False, plot=False)

            r0, r1, r2 = result[0], result[1], result[2]

            # Only skip if PV or Total is bad — allow HV to be NaN
            if np.isnan(r0["pew_meas"]) or np.isnan(r2["pew_meas"]):
                continue

            pew_pv_list.append(r0["pew_meas"])
            vel_pv_list.append(r0["velocity"])
            fwhm_pv_list.append(r0["wmeas"])

            if not np.isnan(r1["pew_meas"]):  # Add HV only if valid
                pew_hv_list.append(r1["pew_meas"])
                vel_hv_list.append(r1["velocity"])
                fwhm_hv_list.append(r1["wmeas"])

            pew_total_list.append(r2["pew_meas"])
        except Exception as e:
            continue

    return {
        "pew_pv_err": np.std(pew_pv_list) if len(pew_pv_list) >= 2 else np.nan,
        "pew_hv_err": np.std(pew_hv_list) if len(pew_hv_list) >= 2 else np.nan,
        "pew_total_err": np.std(pew_total_list) if len(pew_total_list) >= 2 else np.nan,
        "vel_pv_err": np.std(vel_pv_list) if len(vel_pv_list) >= 2 else np.nan,
        "vel_hv_err": np.std(vel_hv_list) if len(vel_hv_list) >= 2 else np.nan,
        "fwhm_pv_err": np.std(fwhm_pv_list) if len(fwhm_pv_list) >= 2 else np.nan,
        "fwhm_hv_err": np.std(fwhm_hv_list) if len(fwhm_hv_list) >= 2 else np.nan,
        "n_valid": len(pew_pv_list)
    }
    
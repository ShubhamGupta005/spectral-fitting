import numpy as np
from scipy.optimize import curve_fit
from copy import deepcopy

from models import mult_gauss_flexible, get_gaussians_paramflex, model_flux
from preprocessing import clean_spectrum, fit_continuum
from plotting import plot_fit
from utils import trapezoid_area, evaluate_hv_snr
from guesses import build_siII_guess, build_cahk_guess



# --- Unified fitting engine ---
def fit_line(spec, redshift, rest1, rest2, window, blue, red, labels, title,
             tie_config, p0, bounds, use_snr_filter=True, plot=True):
    """
    Purpose:
    • Fit a normalized spectrum using a combination of four Gaussians (PV and HV doublet components) plus a constant continuum.
    • Selects between models with and without HV based on Δχ² comparison.
    • Computes velocities, widths (FWHM), and pseudo-equivalent widths (pEW) of PV, HV, and total features.
    • Optionally generates a diagnostic plot of the fit.

    Parameter Descriptions:
    • spec : dict
        Dictionary containing 'wave' and 'flux' arrays of the observed spectrum.
    • redshift : float
        Redshift of the source, used to convert observed wavelengths to rest-frame.
    • rest1 : float
        Rest wavelength of the bluer line in the doublet (e.g., 3934.777 for Ca II K).
    • rest2 : float
        Rest wavelength of the redder line in the doublet (e.g., 3969.591 for Ca II H).
    • window : tuple of float
        Wavelength range (min, max) over which to extract the fitting region.
    • blue, red : float
        Wavelength bounds for the region containing the spectral feature (for masking and normalization).
    • labels : list of str
        Labels for the four Gaussian components (e.g., ["PV K", "PV H", "HV K", "HV H"]).
    • title : str
        Title to be used in the plot.
    • tie_config : dict, optional
        Dictionary of tie settings: {'stied', 'amp_ratio', 'vtied', 'wtied'}.
    • use_snr_filter : bool, optional
        Whether to suppress HV components with low SNR or low amplitude. Default is True.
    • plot : bool, optional
        Whether to show the plot and print fit diagnostics. Default is True.

    Returns:
    • specfeat : list of dict
        List of three dictionaries:
        - PV feature: {"restwl", "velocity", "wmeas" (FWHM), "pew_meas"}
        - HV feature: same keys
        - Total feature: pEW only; velocity and width are set to NaN

    Physics Concept:
    • Spectral absorption features are modeled as Gaussians shifted by Doppler velocity (v/c).
    • Curve fitting minimizes residuals between observed and model flux.
    • Δχ² = χ²_noHV - χ²_full is used to decide whether the HV component improves the fit.
    • pEW is calculated by trapezoidal integration of absorption relative to continuum.
    """
    amp_ratio = tie_config.get("amp_ratio", 1.0)
    vtied = tie_config.get("vtied", True)
    wtied = tie_config.get("wtied", True)

    opz = 1 + redshift
    wave_rest = spec['wave'] / opz
    flux = np.array(spec['flux'])

    mask = (wave_rest >= window[0]) & (wave_rest <= window[1])
    if np.sum(mask) < 10:
        raise RuntimeError("Too few data points in fit window.")

    x = wave_rest[mask]
    y_raw = flux[mask]
    resolution = np.median(np.diff(x))
    win = int(round(20.0 / resolution))
    y_clean = clean_spectrum(x, y_raw, window=win)
    y_norm = y_clean / fit_continuum(x, y_clean, blue, red)
    x_fit = x[(x > blue) & (x < red)]
    y_fit = y_norm[(x > blue) & (x < red)]

    fit_func = lambda x, *params: mult_gauss_flexible(x, *params, rest1=rest1, rest2=rest2,
                                                      amp_ratio=amp_ratio, vtied=vtied, wtied=wtied)
    from copy import deepcopy
    popt_full, _ = curve_fit(fit_func, x_fit, y_fit, p0=p0, bounds=bounds, maxfev=10000)
    chi2_full = np.sum((y_fit - fit_func(x_fit, *popt_full)) ** 2)

    p0_nohv = deepcopy(p0)
    p0_nohv[3] = 0.0
    bounds_nohv = deepcopy(bounds)
    bounds_nohv[0][3] = 0.0
    bounds_nohv[1][3] = 1e-8
    popt_nohv, _ = curve_fit(fit_func, x_fit, y_fit, p0=p0_nohv, bounds=bounds_nohv, maxfev=10000)
    chi2_nohv = np.sum((y_fit - fit_func(x_fit, *popt_nohv)) ** 2)

    delta_chi2 = chi2_nohv - chi2_full
    use_hv = delta_chi2 >= 1.0
    popt = popt_full if use_hv else popt_nohv

    g1, g2, g3, g4 = get_gaussians_paramflex(x_fit, popt, rest1, rest2, amp_ratio, vtied, wtied)
    model_y = model_flux(popt[-1], g1, g2, g3, g4)
    snr_hv = evaluate_hv_snr(x_fit, y_fit, g3, g4)

    vel_pv = popt[1]
    sigma_pv = popt[2]
    amp_hv = popt[3]
    idx = 4
    vel_hv = popt[idx]; idx += 1
    sigma_hv = popt[idx]; idx += 1

    # Compute temporary pEWs before swap check
    pew_pv_temp = trapezoid_area(g1 + g2, x_fit)
    pew_hv_temp = trapezoid_area(g3 + g4, x_fit)
    
    # Check if HV is significant
    hv_present = amp_hv >= 0.05 and pew_hv_temp >= 0.01 and snr_hv >= 2.5
    
    # Only swap if HV is real and PV velocity is more extreme than HV
    if hv_present and abs(vel_pv) > abs(vel_hv):
        print("Swapped PV and HV...")
        # swap everything
        vel_pv, vel_hv = vel_hv, vel_pv
        sigma_pv, sigma_hv = sigma_hv, sigma_pv
        g1, g2, g3, g4 = g3, g4, g1, g2
        pew_pv, pew_hv = pew_hv_temp, pew_pv_temp
        labels = [labels[2], labels[3], labels[0], labels[1]]
    else:
        pew_pv, pew_hv = pew_pv_temp, pew_hv_temp

    pew_total = trapezoid_area(1 - model_y, x_fit)
    if plot:
        print('==========================================================')
        plot_fit(x_fit, y_fit, model_y, g1, g2, g3, g4, labels, title, blue, red)
        print(f"\u0394\u03c7² = {delta_chi2:.2f} — {'HV retained' if use_hv else 'HV suppressed'}")

    specfeat = [
        {"restwl": rest1, "velocity": vel_pv, "wmeas": sigma_pv * 2.3548, "pew_meas": pew_pv},
        {"restwl": rest2, "velocity": vel_hv, "wmeas": sigma_hv * 2.3548, "pew_meas": pew_hv},
        {"restwl": "Total", "velocity": np.nan, "wmeas": np.nan, "pew_meas": pew_total}
    ]

    if amp_hv < 0.05 or pew_hv < 0.01 or (use_snr_filter and snr_hv < 2.5):
        specfeat[1]["velocity"] = np.nan
        specfeat[1]["wmeas"] = np.nan
        specfeat[1]["pew_meas"] = np.nan

    return specfeat



# --- Driver functions ---
def fit_siII(spec, redshift, use_snr_filter=True, plot=True):
    p0, bounds = build_siII_guess()
    return fit_line(spec, redshift, rest1=6347, rest2=6371, window=(6000, 6500), blue=6050, red=6275,
                    labels=["PV 6347", "PV 6371", "HV 6347", "HV 6371"],
                    title="Si II 6355 (6347 + 6371) — PV and HV Components",
                    tie_config={
                        "stied": True, "amp_ratio": 0.564, "vtied": True, "wtied": True
                    },
                    p0=p0, bounds=bounds, use_snr_filter=use_snr_filter, plot=plot)

def fit_cahk(spec, redshift, use_snr_filter=True, plot=True):
    p0, bounds = build_cahk_guess()
    return fit_line(spec, redshift, rest1=3934.777, rest2=3969.591, window=(3550, 4050), blue=3600, red=4000,
                    labels=["PV K", "PV H", "HV K", "HV H"],
                    title="Ca II HK — PV and HV Components",
                    tie_config={
                        "stied": False, "amp_ratio": 1.0, "vtied": True, "wtied": True
                    },
                    p0=p0, bounds=bounds, use_snr_filter=use_snr_filter, plot=plot)



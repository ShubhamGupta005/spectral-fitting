
import numpy as np
def get_gaussians_paramflex(x_fit, p, rest1, rest2, amp_ratio, vtied, wtied):
    """
    Purpose:
    • Construct four Gaussian absorption components representing the PV and HV features of a doublet (e.g., Si II or Ca II HK).
    • Applies velocity shifts and line strength ratios to simulate Doppler-shifted absorption lines for both components.

    Parameter Descriptions:
    • x_fit : np.ndarray
        Wavelength array (rest-frame) over which the model is evaluated.
    • p : list or np.ndarray
        List of fitting parameters: [amp_pv, vel_pv, sigma_pv, amp_hv, vel_hv, sigma_hv, continuum_level].
    • rest1 : float
        Rest wavelength of the bluer component in the doublet.
    • rest2 : float
        Rest wavelength of the redder component in the doublet.
    • amp_ratio : float
        Multiplicative factor between the amplitudes of rest1 and rest2 for both PV and HV (e.g., 0.564 for Si II).
    • vtied : bool
        Whether to tie the velocities of the two lines in each component (always true in this structure).
    • wtied : bool
        Whether to tie the widths of the two lines in each component (always true in this structure).

    Returns:
    • g1 : np.ndarray
        Gaussian profile for PV component at rest1.
    • g2 : np.ndarray
        Gaussian profile for PV component at rest2 (scaled by amp_ratio).
    • g3 : np.ndarray
        Gaussian profile for HV component at rest1.
    • g4 : np.ndarray
        Gaussian profile for HV component at rest2 (scaled by amp_ratio).

    Physics Concept:
    • Uses Doppler shift formula to compute the shifted wavelength λ = λ₀ × (1 + v/c).
    • Assumes optically thin Gaussian absorption for both PV and HV components.
    • Applies a fixed amplitude ratio between doublet lines (based on quantum mechanical oscillator strengths or empirical calibration).
    """
    c = 299792.458
    amp_pv, vel_pv, sigma_pv = p[0], p[1], p[2]
    amp_hv = p[3]
    idx = 4
    vel_hv = p[idx]; idx += 1
    sigma_hv = p[idx]; idx += 1
    cont = p[idx]

    lambda_pv_1 = rest1 * (1 + vel_pv / c)
    lambda_pv_2 = rest2 * (1 + vel_pv / c)
    lambda_hv_1 = rest1 * (1 + vel_hv / c)
    lambda_hv_2 = rest2 * (1 + vel_hv / c)

    g1 = amp_pv * np.exp(-0.5 * ((x_fit - lambda_pv_1) / sigma_pv) ** 2)
    g2 = amp_pv * amp_ratio * np.exp(-0.5 * ((x_fit - lambda_pv_2) / sigma_pv) ** 2)
    g3 = amp_hv * np.exp(-0.5 * ((x_fit - lambda_hv_1) / sigma_hv) ** 2)
    g4 = amp_hv * amp_ratio * np.exp(-0.5 * ((x_fit - lambda_hv_2) / sigma_hv) ** 2)
    return g1, g2, g3, g4

def model_flux(cont, g1, g2, g3, g4):
    """
    Purpose:
    • Construct the full model spectrum by subtracting the combined Gaussian absorption features from a continuum level.
    • Represents how the observed spectrum would appear with PV and HV absorption features applied.

    Returns:
    • model : np.ndarray
        The synthetic model flux: continuum minus the sum of all Gaussian absorption components.

    Physics Concept:
    • Assumes additive Gaussian absorption lines subtracted from the continuum.
    • Models absorption in normalized flux spectra: F_obs = F_cont - ∑(Gaussian_absorptions).
    • This approach reflects the classical picture of line absorption in stellar or supernova atmospheres.
    """
    return cont - (g1 + g2 + g3 + g4)

def mult_gauss_flexible(x, *params, rest1, rest2, amp_ratio=1.0, vtied=True, wtied=True):
    """
    Purpose:
    • Compute the full absorption line model by assembling multiple Gaussians using the provided parameters.
    • This function is used as the model passed into the curve fitting routine (e.g., `scipy.optimize.curve_fit`).
    Returns:
    • model_flux : np.ndarray
        Computed model flux values over `x`, representing the continuum minus the Gaussian absorptions.

    Physics Concept:
    • Uses Doppler shift λ_obs = λ_rest × (1 + v/c) to shift each Gaussian profile based on velocity.
    • Gaussian line shapes are used to simulate thermal or turbulent broadening of absorption features.
    • The absorption model is subtracted from a flat continuum to mimic how real absorption dips appear in normalized spectra.
    """

    g1, g2, g3, g4 = get_gaussians_paramflex(x, params, rest1, rest2, amp_ratio, vtied, wtied)
    return model_flux(params[-1], g1, g2, g3, g4)

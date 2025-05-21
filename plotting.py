
import matplotlib.pyplot as plt
import numpy as np

def plot_fit(x_fit, y_fit, model_y, g1, g2, g3, g4, labels, title, blue, red):
    """
    Plot spectral fitting result including PV and HV Gaussian components.

    Parameters:
    - x_fit : np.ndarray
        Wavelengths in fitting range.
    - y_fit : np.ndarray
        Normalized observed flux.
    - model_y : np.ndarray
        Full model fit to the spectrum.
    - g1 to g4 : np.ndarray
        Individual Gaussian components (PV/HV).
    - labels : list of str
        Labels for Gaussian components.
    - title : str
        Plot title.
    - blue, red : float
        Fitting window limits.
    """
    plt.figure(figsize=(10, 5))
    plt.plot(x_fit, y_fit, color='gray',label='Normalized Flux')
    plt.plot(x_fit, model_y, 'g-', label='Total Fit')
    plt.plot(x_fit, 1 - g1, 'r--', label=labels[0])
    plt.plot(x_fit, 1 - g2, 'r:', label=labels[1])
    plt.plot(x_fit, 1 - g3, 'b--', label=labels[2])
    plt.plot(x_fit, 1 - g4, 'b:', label=labels[3])
    plt.fill_between(x_fit, 1, 1 - g1 - g2, where=(g1 + g2) > 0, color='red', alpha=0.2, label='pEW PV')
    plt.fill_between(x_fit, 1, 1 - g3 - g4, where=(g3 + g4) > 0, color='blue', alpha=0.2, label='pEW HV')
    plt.xlabel('Rest-frame Wavelength (Å)')
    plt.ylabel('Normalized Flux')
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.xlim(blue, red)
    plt.tight_layout()
    plt.show()
# demo.py — Simple CLI demo for spectral fitting

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from fitting import fit_siII  # Or use fit_cahk for Ca II

# Load example spectrum
import sys

csv_path = sys.argv[1] if len(sys.argv) > 1 else 'example_spectrum.csv'
wave = df['Wavelength'].values
flux = df['Flux'].values
spec = {'wave': wave, 'flux': flux}
z = 0.05  # Example redshift

# Run the fit
result = fit_siII(spec, redshift=z, plot=True)

# Print output
print("\nFit Results:")
for label, r in zip(['PV', 'HV', 'Total'], result):
    vel = r['velocity'] if not np.isnan(r['velocity']) else '—'
    fwhm = r['wmeas'] if not np.isnan(r['wmeas']) else '—'
    print(f"{label}: pEW = {r['pew_meas']:.2f} Å, Velocity = {vel} km/s, FWHM = {fwhm} km/s")

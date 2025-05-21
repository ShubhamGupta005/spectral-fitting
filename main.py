import pandas as pd
import os

# --- CONFIGURATION ---
# metadata_file = '/Users/gupta/Downloads/Project/DES/DES_DATA/02-DATA_SPECTRA/data_parameters/updated_para.csv'
metadata_file = '/Users/gupta/Downloads/Project/DES/final_txt/DES_para_ca.txt'
data_dir = '/Users/gupta/Downloads/Project/DES/DES_DATA/02-DATA_SPECTRA/DES-SN3YR_SNIa_corrected/DES_smooth'

line_type = "cahk"  # or "CaHK"
monte_carlo = True
n_max = 5 # set to an integer to limit number of rows

# --- LOAD METADATA ---
meta = pd.read_csv(metadata_file)
if n_max:
    meta = meta.head(n_max)

# --- PROCESS EACH SPECTRUM ---
for _, row in meta.iterrows():
    try:
        filename, snid, z = row['filename'], row['SNID'], row['REDSHIFT']
        filepath = os.path.join(data_dir, filename)
        df = pd.read_csv(filepath)
        wave, flux = df['Wavelength'].values, df['Flux'].values
        spec = {'wave': wave, 'flux': flux}

        # Run main fit
        if line_type.lower() == 'si':
            result = fit_siII(spec, redshift=z)
        elif line_type.lower() == 'cahk':
            result = fit_cahk(spec, redshift=z)
        else:
            raise ValueError("Invalid line_type: choose 'si' or 'CaHK'.")

        # Run optional MC error analysis
        if monte_carlo:
            mc_errs = run_monte_carlo(spec, redshift=z, line_type=line_type)

        # Print output
        print(f"SNID: {snid}")
        labels = ['PV', 'HV', 'Total']
        for r, label in zip(result, labels):
            if label == "Total":
                print(f"  {label}: pEW = {r['pew_meas']:.2f} Å")
            else:
                print(f"  {label}: pEW = {r['pew_meas']:.2f} Å, Vel_{label} = {r['velocity'] if not np.isnan(r['velocity']) else '—'} km/s, FWHM = {r['wmeas'] if not np.isnan(r['wmeas']) else '—'} km/s")

        if monte_carlo and mc_errs["n_valid"] >= 2:
            print("  Monte Carlo Errors:")
            print(f"    PV: pEW ±{mc_errs['pew_pv_err']:.2f}, Vel ±{mc_errs['vel_pv_err']:.1f}, FWHM ±{mc_errs['fwhm_pv_err']:.2f}")
            print(f"    HV: pEW ±{mc_errs['pew_hv_err']:.2f}, Vel ±{mc_errs['vel_hv_err']:.1f}, FWHM ±{mc_errs['fwhm_hv_err']:.2f} (from {mc_errs['n_valid']} trials)")
            print(f"    Total: pEW ±{mc_errs['pew_total_err']:.2f}")
        elif monte_carlo:
            print(f"  Monte Carlo skipped: only {mc_errs['n_valid']} valid trial(s)")
        print()

    except Exception as e:
        print(f"Failed for {filename} ({snid}): {e}")
# Spectral Fitting Package for Type Ia Supernovae

This Python package fits high-velocity (HV) and photospheric-velocity (PV) components in Type Ia supernova spectra, focusing on Si II 6355 and Ca II HK doublets. It supports multi-Gaussian absorption modeling, pEW measurements, and Monte Carlo error estimation.

---

## 📁 Project Structure

```
spectral_fitting/
├── fitting.py          # Core fitting logic and Gaussian fitting engine
├── models.py           # Gaussian profile construction (PV, HV components)
├── preprocessing.py    # Spectrum cleaning and continuum fitting
├── plotting.py         # Fit visualization (flux, models, pEW shading)
├── utils.py            # SNR and trapezoidal integration
├── guesses.py          # Initial guess builders for Si II and Ca II
├── montecarlo.py       # Monte Carlo error estimation engine
├── main.py             # Batch fitting driver using metadata
├── requirements.txt    # Package dependencies
└── README.md           # This file
```

---

## ⚙️ Installation

1. Clone the repo or download the files:

```bash
git clone https://github.com/your_username/spectral-fitting.git
cd spectral-fitting
```

2. Install dependencies:

```bash
pip install -r requirements.txt
```

---

## 🚀 Usage

### Fit a single spectrum:
```python
from fitting import fit_siII

spec = {
    'wave': your_wavelength_array,
    'flux': your_flux_array
}
result = fit_siII(spec, redshift=0.05)
```

### Run batch fitting:
Edit `main.py` to set:
- `metadata_file`: CSV/TXT with columns `filename`, `SNID`, `REDSHIFT`
- `data_dir`: folder containing the spectral files
- `line_type`: `'si'` or `'cahk'`

Then run:
```bash
python main.py
```

---

## 📊 Output

For each spectrum:
- pEW, velocity, and FWHM for PV, HV, and total components
- Optional Monte Carlo errors if enabled
- Diagnostic fit plot showing absorption features

---

## 📦 Dependencies

- numpy
- scipy
- pandas
- matplotlib

Install with:

```bash
pip install -r requirements.txt
```

---

## 🧪 Line Types Supported

| Line     | Rest Wavelengths (Å) | Notes                          |
|----------|----------------------|--------------------------------|
| Si II    | 6347, 6371           | Amp ratio = 0.564, vtied = True |
| Ca II HK | 3934.777, 3969.591   | Independent lines              |

---

## 📄 License

MIT License

---

## 🤝 Credits

Developed by Ray for Type Ia Supernova Spectral Analysis.

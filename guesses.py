# --- Initial guess builder functions for each feature ---
def build_siII_guess():
    vel_pv = -12000
    vel_hv = vel_pv - 6000
    vel_hv_min = vel_pv - 4500

    p0 = [0.05, vel_pv, 10.0, 0.02, max(vel_hv, vel_hv_min), 10.0, 1.0]
    bounds_low = [0, -20000, 2.0, 0, vel_hv_min, 2.0, 0.98]
    bounds_high = [1, -5000, 20.0, 1, -5000, 20.0, 1.02]
    return p0, (bounds_low, bounds_high)

def build_cahk_guess():
    p0 = [0.05, -12000, 10.0, 0.02, -18000, 10.0, 1.0]
    bounds_low = [0, -20000, 2.0, 0, -20000, 2.0, 0.98]
    bounds_high = [1, -5000, 20.0, 1, -5000, 20.0, 1.02]
    return p0, (bounds_low, bounds_high)

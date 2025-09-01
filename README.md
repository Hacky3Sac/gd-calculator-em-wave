
[pectral_line_id_tool.py](https://github.com/user-attachments/files/22084130/pectral_line_id_tool.py)
# Full Spectral Line ID from a single wavelength (vacuum)
# Edit LAMBDA_NM and run. No external libraries required.

# ---- EDIT THIS LINE ----
LAMBDA_NM = 254.02137   # example: Hydrogen Hα (nm)
# ------------------------

import math

# Physical constants (SI)
c   = 299_792_458.0                # m/s
h   = 6.62607015e-34               # J·s
kB  = 1.380649e-23                 # J/K
eC  = 1.602176634e-19              # C  (J <-> eV)
G   = 6.67430e-11                  # m^3·kg^-1·s^-2 (gravitational constant, 2018 CODATA)

# Convert input
lam = LAMBDA_NM * 1e-9             # wavelength in meters

# Core wave/photon relations
f     = c / lam                    # Hz
omega = 2 * math.pi * f            # rad/s
C_circ= lam                        # m  (C ≡ λ)
r     = lam / (2 * math.pi)        # m
k_num = 2 * math.pi / lam          # 1/m
T_per = 1.0 / f                    # s
E_J   = h * f                      # J
E_eV  = E_J / eC                   # eV
p     = h / lam                    # kg·m/s
m_eq  = E_J / (c**2)               # kg  (energy-equivalent, not photon rest mass)
v     = c                          # m/s  (vacuum EM wave)
F_Er  = E_J / r                    # N   (J/m)
F_g   = G * (m_eq**2) / (r**2)     # N   (using m_eq; illustrative only)
T_eff = E_J / kB                   # K   (energy scale)

# Pretty print
def pr(name, val, unit):
    print(f"{name:<35} {val: .6e} {unit}")

print(f"\nFull Spectral Line ID (λ = {LAMBDA_NM:.6f} nm, vacuum)\n")
pr("Speed of light (c)",            c,      "[m/s]")
pr("Wavelength (λ)",                lam,    "[m]")
pr("Frequency (f = c/λ)",           f,      "[Hz]")
pr("Angular velocity (ω = 2πf)",    omega,  "[rad/s]")
pr("Circumference (C = λ)",         C_circ, "[m]")
pr("Radius (r = λ/2π)",             r,      "[m]")
pr("Wave number (k = 2π/λ)",        k_num,  "[m^-1]")
pr("Period (T = 1/f)",              T_per,  "[s]")
pr("Energy (E = hf)",               E_J,    "[J]")
pr("Energy (eV)",                   E_eV,   "[eV]")
pr("Momentum (p = h/λ)",            p,      "[kg*m/s]")
pr("Equivalent mass (m = E/c^2)",    m_eq,   "[kg]")
pr("Velocity (v)",                  v,      "[m/s]")
pr("Force eq. (E/r)",               F_Er,   "[N]")
pr("Grav. force at r (G m^2/r^2)",    F_g,    "[N]")
pr("Temperature eq. (E/kB)",        T_eff,  "[K]")



# Optional quick checks (relative errors; should be ~0)
rel = lambda a,b: 0.0 if a==0 else abs(a-b)/abs(a)
print("\nSanity checks (relative error):")
print(f"  ω vs 2πf       : {rel(omega, 2*math.pi*f):.3e}")
print(f"  ω vs c/r       : {rel(omega, c/r):.3e}")
print(f"  E vs h f       : {rel(E_J, h*f):.3e}")
print(f"  E vs p c       : {rel(E_J, p*c):.3e}\n")

# Notes:
# - 'Equivalent mass' uses E/c^2; photons have zero rest mass.
# - 'Force eq.' E/r and 'Grav. force' G m^2 / r^2 are dimensional constructs,
#   not forces acting on a specific body; use with physical care.

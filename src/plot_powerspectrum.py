import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from astropy import constants as const
from astropy import units as u
from scipy.signal import argrelextrema as extrema

plt.style.use("bmh")

fonts = {
    "font.family": "serif",
    "axes.labelsize": 18,
    "font.size": 12,
    "legend.fontsize": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
}

plt.rcParams.update(fonts)

# Loading background cosmology data
cosmo_data = np.loadtxt("cosmology.txt")
x_cosmo     = cosmo_data[:, 0]
eta_cosmo = (cosmo_data[:, 1] * u.m).to(u.Mpc)
Hp_cosmo    = (cosmo_data[:, 2] / u.s).to(u.km / (u.s * u.Mpc))
OmegaB      = cosmo_data[:, 3]
OmegaCDM    = cosmo_data[:, 4]
OmegaLambda = cosmo_data[:, 5]
OmegaR      = cosmo_data[:, 6]
Omega_sum   = OmegaB + OmegaCDM + OmegaLambda + OmegaR
Omega_m     = OmegaCDM + OmegaB

eq_index  = np.argmin(np.abs(OmegaR - Omega_m))
k_horizon = 1 / eta_cosmo 
k_eq      = k_horizon[eq_index]

# Loading power spectrum data
data = np.loadtxt("cells.txt")
data_test = np.loadtxt("cells_test.txt")
ells = data[:, 0]
Cells = data[:, 1]

# Cells for different matter-energy content 
data_test2 = np.loadtxt("cells_test2.txt")
Cells_test = data_test[:, 1]
Cells_test2 = data_test2[:, 1]

data1 = np.loadtxt("matter_power_spec.txt")
k = data1[:, 0] * 1 / u.Mpc
P = data1[:, 1] * u.Mpc ** 3

# Matter power spectrum for different matter-energy content
data1_test = np.loadtxt("matter_power_spec_test.txt")
P_test = data1_test[:, 1] * u.Mpc ** 3

Theta_ell_peak_1 = data1[:, 2]
Theta_ell_peak_2 = data1[:, 3]
Theta_ell_peak_3 = data1[:, 4]
Theta_ell_peak_4 = data1[:, 5]

Theta_ell_sq_peak_1 = data1[:, 6] * u.Mpc
Theta_ell_sq_peak_2 = data1[:, 7] * u.Mpc
Theta_ell_sq_peak_3 = data1[:, 8] * u.Mpc
Theta_ell_sq_peak_4 = data1[:, 9] * u.Mpc

Cell_peaks = extrema(Cells, np.greater)[0]
Cell_troughs = extrema(Cells, np.less)[0]


print(f"k_* = {Cell_peaks[0] / eta_cosmo[-1].value} h/Mpc")

# Plotting temperature power spectrum with peaks and troughs marking
fig, ax = plt.subplots(3, 1, figsize=[1.5 * 7.1014, 1.5 * 7.1014 / 1.618])
ax[0].axvline(ells[Cell_peaks[0]], color = "r", label = rf"$\ell = {ells[Cell_peaks[0]]}$")
ax[0].axvline(ells[Cell_peaks[1]], color = "g", label = rf"$\ell = {ells[Cell_peaks[1]]}$")
ax[0].axvline(ells[Cell_peaks[2]], color = "b", label = rf"$\ell = {ells[Cell_peaks[2]]}$")
ax[0].axvline(ells[Cell_peaks[3]], color = "m", label = rf"$\ell = {ells[Cell_peaks[3]]}$")
ax[0].axvline(ells[Cell_troughs[0]], color = "r", linestyle = ":", label = rf"$\ell = {ells[Cell_troughs[0]]}$")
ax[0].axvline(ells[Cell_troughs[1]], color = "g", linestyle = ":", label = rf"$\ell = {ells[Cell_troughs[1]]}$")
ax[0].axvline(ells[Cell_troughs[2]], color = "b", linestyle = ":", label = rf"$\ell = {ells[Cell_troughs[2]]}$")
ax[0].axvline(ells[Cell_troughs[3]], color = "m", linestyle = ":", label = rf"$\ell = {ells[Cell_troughs[3]]}$")
ax[0].axvline(ells[Cell_troughs[4]], color = "orange", linestyle = "--", label = rf"$\ell = {ells[Cell_troughs[4]]}$")
ax[0].plot(ells, Cells)
ax[0].set_yscale("log")
ax[0].set_xscale("log")
#ax[0].set_xlabel(r"$\ell$")
ax[0].set_ylabel(r"$\frac{\ell(\ell + 1)}{2\pi}C_\ell$ $[\mu \mathrm{K}^2]$")
ax[0].legend(loc = 0, fontsize = 9.5, ncol = 4)

# Plotting tenmperature power spectrum for different parameter combos
ax[1].plot(ells, Cells_test, linestyle = "--", label = r"$(\Omega_{CDM0}, \Omega_{B0}) = (0.173, 0.1)$")
ax[1].plot(ells, Cells_test2, linestyle = "--", label = r"$(\Omega_{CDM0}, \Omega_{B0}) = (0.049, 0.224)$")
ax[1].plot(ells, Cells, label = r"$(\Omega_{CDM0}, \Omega_{B0}) = (0.224, 0.049)$")
ax[1].set_yscale("log")
ax[1].set_xscale("log")
ax[1].set_xlabel(r"$\ell$")
ax[1].set_ylabel(r"$\frac{\ell(\ell + 1)}{2\pi}C_\ell$ $[\mu \mathrm{K}^2]$")
ax[1].legend(loc = 0, fontsize = 9.5)

# Plotting matter power spectrum
ax[2].plot(k, P , label = r"$P_M(k), (\Omega_{CDM0}, \Omega_{B0}) = (0.224, 0.049)$")
ax[2].plot(k, P_test, color = "cornflowerblue", linestyle = "--", label = r"$P_M(k), (\Omega_{CDM0}, \Omega_{B0}) = (0.173, 0.1)$")
ax[2].axvline(k_eq.value, color = "r", linestyle = "--", label = rf"$k_\mathrm{{eq}} = {k_eq.value:.2g}\mathrm{{h/Mpc}}$")
ax[2].set_yscale("log")
ax[2].set_xscale("log")
ax[2].set_xlabel(r"$k$ $[h/\mathrm{Mpc}^{-1}]$")
ax[2].set_ylabel(r"$P_M(k)$ $[(\mathrm{Mpc}/h)^3]$")
ax[2].legend(loc = 0)
fig.tight_layout(pad = 0.1)
plt.savefig("../doc/Figures/Cell.pdf")


fig1, ax1 = plt.subplots(2, 1, figsize=[1.5 * 7.1014, 1.5 * 7.1014 / 1.618])
# Plotting integrands for different scales ell
ax1[0].plot(k, Theta_ell_peak_1, label=rf"$\ell = 10$", color = "m")
ax1[0].plot(k, Theta_ell_peak_2, label=rf"$\ell = 25$", color = "b")
ax1[0].plot(k, Theta_ell_peak_3, label=rf"$\ell = 50$", color = "r")
ax1[0].plot(k, Theta_ell_peak_4, label=rf"$\ell = 100$", color = "g")
ax1[0].set_ylabel(r"$\Theta_\ell(k)$")
ax1[0].set_xlabel(r"$k$ $[h/\mathrm{Mpc}^{-1}]$")
ax1[0].legend(loc = 0)
ax1[0].set_xlim(0, 0.05)
ax1[0].set_xscale("linear")
ax1[0].set_yscale("linear")

# Plotting integrands squared for different scales ell
ax1[1].plot(k, Theta_ell_sq_peak_1, label=rf"$\ell = 10$", color = "m")
ax1[1].plot(k, Theta_ell_sq_peak_2, label=rf"$\ell = 25$", color = "b")
ax1[1].plot(k, Theta_ell_sq_peak_3, label=rf"$\ell = 50$", color = "r")
ax1[1].plot(k, Theta_ell_sq_peak_4, label=rf"$\ell = 100$", color = "g")
ax1[1].set_ylabel(r"$\Theta_\ell^2(k)/k$ $[Mpc / h]$")
ax1[1].set_xlabel(r"$k$ $[h/\mathrm{Mpc}^{-1}]$")
ax1[1].legend(loc = 0)
ax1[1].set_xlim(0, 0.05)
ax1[1].set_xscale("linear")
ax1[1].set_yscale("linear")
fig1.tight_layout()
plt.savefig("../doc/Figures/Transfer_func.pdf")


plt.show()

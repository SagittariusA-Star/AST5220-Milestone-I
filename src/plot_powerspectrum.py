import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from astropy import constants as const
from astropy import units as u

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

data = np.loadtxt("cells.txt")
ells = data[:, 0]
Cells = data[:, 1]

data1 = np.loadtxt("matter_power_spec.txt")
k = data1[:, 0] * 1 / u.Mpc
P = data1[:, 1] * u.Mpc ** 3
Theta_ell_peak_1 = data1[:, 2]
Theta_ell_peak_2 = data1[:, 3]
Theta_ell_peak_3 = data1[:, 4]
Theta_ell_peak_4 = data1[:, 5]

Theta_ell_trough_1 = data1[:, 5]
Theta_ell_trough_2 = data1[:, 6]
Theta_ell_trough_3 = data1[:, 7]
Theta_ell_trough_4 = data1[:, 8]


Theta_ell_sq_peak_1 = data1[:, 9] * u.Mpc
Theta_ell_sq_peak_2 = data1[:, 10] * u.Mpc
Theta_ell_sq_peak_3 = data1[:, 11] * u.Mpc
Theta_ell_sq_peak_4 = data1[:, 12] * u.Mpc

Theta_ell_sq_trough_1 = data1[:, 13] * u.Mpc
Theta_ell_sq_trough_2 = data1[:, 14] * u.Mpc
Theta_ell_sq_trough_3 = data1[:, 15] * u.Mpc
Theta_ell_sq_trough_4 = data1[:, 16] * u.Mpc

fig, ax = plt.subplots(2, 1, figsize=[1.5 * 7.1014, 1.5 * 7.1014 / 1.618])
ax[0].plot(ells, Cells)
ax[0].set_yscale("log")
ax[0].set_xscale("log")
ax[0].set_xlabel(r"$\ell$")
ax[0].set_ylabel(r"$\frac{\ell(\ell + 1)}{2\pi}C_\ell$ $[\mu \mathrm{K}^2]$")

ax[1].plot(k, P , label = r"$P_M(k)$")
ax[1].axvline(k_eq.value, color = "r", linestyle = "--", label = rf"$k_\mathrm{{eq}} = {k_eq.value:.2g}\mathrm{{h/Mpc}}$")
ax[1].set_yscale("log")
ax[1].set_xscale("log")
ax[1].set_xlabel(r"$k$ $[h/\mathrm{Mpc}]^{-1}$")
ax[1].set_ylabel(r"$P_M(k)$ $[(\mathrm{Mpc}/h)^3]$")
ax[1].legend(loc = 0)
fig.tight_layout()
plt.savefig("../doc/Figures/Cell.pdf")


fig1, ax1 = plt.subplots(2, 2, figsize=[1.5 * 7.1014, 1.5 * 7.1014 / 1.618])

# Monopole moments
ax1[0, 0].plot(k, Theta_ell_peak_1, label=rf"$\ell = 200$", color = "m")
ax1[0, 0].plot(k, Theta_ell_peak_2, label=rf"$\ell = 480$", color = "b")
ax1[0, 0].plot(k, Theta_ell_peak_3, label=rf"$\ell = 725$", color = "r")
ax1[0, 0].plot(k, Theta_ell_peak_4, label=rf"$\ell = 1000$", color = "g")
ax1[0, 0].set_ylabel(r"$\Theta_\ell(k)$")
ax1[0, 0].set_xlabel(r"$k [h/\mathrm{Mpc}^{-1}]$")
ax1[0, 0].legend(loc = 0)
ax1[0, 0].set_xscale("linear")
ax1[0, 0].set_yscale("linear")

ax1[0, 1].plot(k, Theta_ell_trough_1, label=rf"$\ell = 8$", color = "m")
ax1[0, 1].plot(k, Theta_ell_trough_2, label=rf"$\ell = 370$", color = "b")
ax1[0, 1].plot(k, Theta_ell_trough_3, label=rf"$\ell = 590$", color = "r")
ax1[0, 1].plot(k, Theta_ell_trough_4, label=rf"$\ell = 896$", color = "g")
ax1[0, 1].set_ylabel(r"$\Theta_\ell(k)$")
ax1[0, 1].set_xlabel(r"$k [h/\mathrm{Mpc}^{-1}]$")
ax1[0, 1].legend(loc = 0)
ax1[0, 1].set_xscale("linear")
ax1[0, 1].set_yscale("linear")

ax1[1, 0].plot(k, Theta_ell_sq_peak_1, label=rf"$\ell = 200$", color = "m")
ax1[1, 0].plot(k, Theta_ell_sq_peak_2, label=rf"$\ell = 480$", color = "b")
ax1[1, 0].plot(k, Theta_ell_sq_peak_3, label=rf"$\ell = 725$", color = "r")
ax1[1, 0].plot(k, Theta_ell_sq_peak_4, label=rf"$\ell = 1000$", color = "g")
ax1[1, 0].set_ylabel(r"$\Theta_\ell^2(k)/k$ $[Mpc / h]")
ax1[1, 0].set_xlabel(r"$k [h/\mathrm{Mpc}^{-1}]$")
ax1[1, 0].legend(loc = 0)
ax1[1, 0].set_xscale("linear")
ax1[1, 0].set_yscale("linear")

ax1[1, 1].plot(k, Theta_ell_sq_trough_1, label=rf"$\ell = 8$", color = "m")
ax1[1, 1].plot(k, Theta_ell_sq_trough_2, label=rf"$\ell = 370$", color = "b")
ax1[1, 1].plot(k, Theta_ell_sq_trough_3, label=rf"$\ell = 590$", color = "r")
ax1[1, 1].plot(k, Theta_ell_sq_trough_4, label=rf"$\ell = 896$", color = "g")
ax1[1, 1].set_ylabel(r"$\Theta_\ell^2(k)/k$ $[Mpc / h]")
ax1[1, 1].set_xlabel(r"$k [h/\mathrm{Mpc}^{-1}]$")
ax1[1, 1].legend(loc = 0)
ax1[1, 1].set_xscale("linear")
ax1[1, 1].set_yscale("linear")

fig1.tight_layout()
plt.savefig("../doc/Figures/Transfer_func.pdf")


plt.show()

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


data = np.loadtxt("cells.txt")
ells = data[:, 0]
Cells = data[:, 1]

data1 = np.loadtxt("matter_power_spec.txt")
k = (data1[:, 0] / u.m).to(1 / u.Mpc)
P = (data1[:, 1] * u.m ** 3).to(u.Mpc ** 3)
Theta_ell_peak_1 = data1[:, 2]
Theta_ell_peak_2 = data1[:, 3]
Theta_ell_peak_3 = data1[:, 4]
Theta_ell_peak_4 = data1[:, 5]

Theta_ell_trough_1 = data1[:, 5]
Theta_ell_trough_2 = data1[:, 6]
Theta_ell_trough_3 = data1[:, 7]
Theta_ell_trough_4 = data1[:, 8]


Theta_ell_sq_peak_1 = data1[:, 9]
Theta_ell_sq_peak_2 = data1[:, 10]
Theta_ell_sq_peak_3 = data1[:, 11]
Theta_ell_sq_peak_4 = data1[:, 12]

Theta_ell_sq_trough_1 = data1[:, 13]
Theta_ell_sq_trough_2 = data1[:, 14]
Theta_ell_sq_trough_3 = data1[:, 15]
Theta_ell_sq_trough_4 = data1[:, 16]

fig = plt.figure(figsize=[1.5 * 7.1014, 1.5 * 7.1014 / 1.618])
ax = plt.subplot(111)
ax.plot(ells, Cells)
ax.set_yscale("log")
ax.set_xscale("log")
plt.savefig("../doc/Figures/Cell.pdf")


fig1, ax1 = plt.subplots(2, 2, figsize=[1.5 * 7.1014, 1.5 * 7.1014 / 1.618])

# Monopole moments
ax1[0, 0].plot(k, Theta_ell_peak_1, label=rf"$\ell = 200$", color = "m")
ax1[0, 0].plot(k, Theta_ell_peak_2, label=rf"$\ell = 480$", color = "b")
ax1[0, 0].plot(k, Theta_ell_peak_3, label=rf"$\ell = 725$", color = "r")
ax1[0, 0].plot(k, Theta_ell_peak_4, label=rf"$\ell = 1000$", color = "g")
ax1[0, 0].set_ylabel(r"$\Theta_\ell(k)$")
ax1[0, 0].set_xlabel(r"$k [\mathrm{Mpc}^{-1}]$")
ax1[0, 0].legend(loc = 0)
ax1[0, 0].set_xscale("log")
ax1[0, 0].set_yscale("linear")

ax1[0, 1].plot(k, Theta_ell_trough_1, label=rf"$\ell = 8$", color = "m")
ax1[0, 1].plot(k, Theta_ell_trough_2, label=rf"$\ell = 370$", color = "b")
ax1[0, 1].plot(k, Theta_ell_trough_3, label=rf"$\ell = 590$", color = "r")
ax1[0, 1].plot(k, Theta_ell_trough_4, label=rf"$\ell = 896$", color = "g")
ax1[0, 1].set_ylabel(r"$\Theta_\ell(k)$")
ax1[0, 1].set_xlabel(r"$k [\mathrm{Mpc}^{-1}]$")
ax1[0, 1].legend(loc = 0)
ax1[0, 1].set_xscale("log")
ax1[0, 1].set_yscale("linear")

ax1[1, 0].plot(k, Theta_ell_sq_peak_1, label=rf"$\ell = 200$", color = "m")
ax1[1, 0].plot(k, Theta_ell_sq_peak_2, label=rf"$\ell = 480$", color = "b")
ax1[1, 0].plot(k, Theta_ell_sq_peak_3, label=rf"$\ell = 725$", color = "r")
ax1[1, 0].plot(k, Theta_ell_sq_peak_4, label=rf"$\ell = 1000$", color = "g")
ax1[1, 0].set_ylabel(r"$\Theta_\ell^2(k)/k$")
ax1[1, 0].set_xlabel(r"$k [\mathrm{Mpc}^{-1}]$")
ax1[1, 0].legend(loc = 0)
ax1[1, 0].set_xscale("log")
ax1[1, 0].set_yscale("linear")

ax1[1, 1].plot(k, Theta_ell_sq_trough_1, label=rf"$\ell = 8$", color = "m")
ax1[1, 1].plot(k, Theta_ell_sq_trough_2, label=rf"$\ell = 370$", color = "b")
ax1[1, 1].plot(k, Theta_ell_sq_trough_3, label=rf"$\ell = 590$", color = "r")
ax1[1, 1].plot(k, Theta_ell_sq_trough_4, label=rf"$\ell = 896$", color = "g")
ax1[1, 1].set_ylabel(r"$\Theta_\ell^2(k)/k$")
ax1[1, 1].set_xlabel(r"$k [\mathrm{Mpc}^{-1}]$")
ax1[1, 1].legend(loc = 0)
ax1[1, 1].set_xscale("log")
ax1[1, 1].set_yscale("linear")

fig1.tight_layout()

fig3 = plt.figure(figsize=[1.5 * 7.1014, 1.5 * 7.1014 / 1.618])
ax3 = plt.subplot(111)
ax3.plot(k, P)
ax3.set_yscale("log")
ax3.set_xscale("log")
ax3.set_xlabel("k")
ax3.set_ylabel("P_M")

plt.show()

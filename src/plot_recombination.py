import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from astropy import constants as const
from astropy import units as u

plt.style.use("bmh")

fonts = {
    "font.family": "serif",
    "axes.labelsize": 12,
    "font.size": 10,
    "legend.fontsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
}

plt.rcParams.update(fonts)

# Loading data from file
data = np.loadtxt("recombination.txt")
x = data[:, 0]
Xe = data[:, 1]
ne = data[:, 2]
tau = data[:, 3]
dtaudx = data[:, 4]
ddtaudxdx = data[:, 5]
g_tilde = data[:, 6]
dg_tildedx = data[:, 7]
ddg_tildeddx = data[:, 8]

fig, ax = plt.subplots(2, 2, figsize=[7.1014, 7.1014 / 1.618])
ax[0, 0].plot(x, Xe, label=r"$X_e(x)$")
ax[0, 0].legend()
ax[0, 0].set_xlabel(r"$x = \log (a)$")
ax[0, 0].set_ylabel(r"$X_e \approx n_e / n_H$")
ax[0, 0].set_yscale("log")
#ax[0, 0].set_xlim(x[0], x[np.where[x >= 0]])

ax[0, 1].plot(x, ne, label=r"$n_e(x)$")
ax[0, 1].legend()
ax[0, 1].set_xlabel(r"$x = \log (a)$")
ax[0, 1].set_ylabel(r"$n_e [\mathrm{m^{-3}}]$")
ax[0, 1].set_yscale("log")
#ax[0, 1].set_xlim(x[0], x[np.where[x >= 0]])

ax[1, 0].plot(x, tau, label=r"$\tau(x)$")
ax[1, 0].plot(x, - dtaudx, label=r"$-\tau'(x)$", linestyle = "--")
ax[1, 0].plot(x, ddtaudxdx, label=r"$\tau''(x)$", linestyle = "-.")
ax[1, 0].legend()
ax[1, 0].set_xlabel(r"$x = \log (a)$")
ax[1, 0].set_ylim(1e-8, 1e8)
ax[1, 0].set_xlim(-12, 0)
ax[1, 0].set_yscale("log")

ax[1, 1].plot(x, g_tilde, label=r"$\tilde{g}(x)$")
ax[1, 1].set_xlim(-12, 0)


ax[1, 1].legend()
ax[1, 1].set_xlabel(r"$x = \log (a)$")
#ax[1, 1].set_xlim(x[0], x[np.where[x >= 0]])

fig1, ax1 = plt.subplots(1, 2, figsize=[7.1014, 7.1014 / 1.618])
ax1[0].plot(x, dg_tildedx, label=r"$\tilde{g}'(x)$")
ax1[0].set_xlabel(r"$x = \log (a)$")
ax1[0].legend()

ax1[1].plot(x, ddg_tildeddx, label=r"$\tilde{g}''(x)$")
ax1[1].legend()
ax1[1].set_xlabel(r"$x = \log (a)$")

plt.show()


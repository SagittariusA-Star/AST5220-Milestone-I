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
    "legend.fontsize": 18,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
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

# Computing printout data
g_integral = np.trapz(g_tilde, x = x)
log_rel_error = np.log10(np.abs(g_integral - 1))
x_rec = x[np.where(g_tilde == g_tilde.max())][0]
a_rec = np.exp(x_rec)
z_rec = 1 / a_rec - 1
g_max = g_tilde.max()
tau_transparent = tau[np.abs(tau - 1).argmin()]
x_transparent = x[np.abs(tau - 1).argmin()]

# Printout of interesting information
print("----------------------------Some interesting quantities------------------------")
print("Integral of g_tilde: {0}, log rel error: {1}".format(g_integral, log_rel_error))
print("Maximum of g_tilde: x = {0}, g_tilde = {1}".format(x_rec, g_max))
print("Optical depth: tau = {0} at x = {1}".format(tau_transparent, x_transparent))
print("Redshift at recombination: z = {0}".format(z_rec))
print("-------------------------------------------------------------------------------")

# Generating plots
fig, ax = plt.subplots(2, 2, figsize=[1.5 * 7.1014, 1.5 * 7.1014 / 1.618])

# Plotting electron fraction
ax[0, 0].plot(x, Xe, label=r"$X_e(x)$")
ax[0, 0].legend()
ax[0, 0].set_xlabel(r"$x = \log (a)$")
ax[0, 0].set_ylabel(r"$X_e \approx n_e / n_H$")
ax[0, 0].set_yscale("log")
ax[0, 0].set_xlim(-12, 0)

# Plotting electron density
ax[0, 1].plot(x, ne, label=r"$n_e(x)$")
ax[0, 1].legend()
ax[0, 1].set_xlabel(r"$x = \log (a)$")
ax[0, 1].set_ylabel(r"$n_e [\mathrm{m^{-3}}]$")
ax[0, 1].set_yscale("log")
ax[0, 1].set_xlim(-12, 0)

# Plotting optical depth
ax[1, 0].plot(x, tau, label=r"$\tau(x)$")
ax[1, 0].plot(x, - dtaudx, label=r"$-\tau'(x)$", linestyle = "--")
ax[1, 0].plot(x, ddtaudxdx, label=r"$\tau''(x)$", linestyle = "-.")
ax[1, 0].legend()
ax[1, 0].set_xlabel(r"$x = \log (a)$")
ax[1, 1].set_ylabel(r"$\tau(x)$")
ax[1, 0].set_ylim(1e-8, 1e8)
ax[1, 0].set_xlim(-12, 0)
ax[1, 0].set_yscale("log")

# Plotting visibility function
ax[1, 1].plot(x, g_tilde, label=r"$\tilde{g}(x)$")
ax[1, 1].set_xlim(-12, 0)
ax[1, 1].legend()
ax[1, 1].set_ylabel(r"$\tilde{g}(x)$")
ax[1, 1].set_xlabel(r"$x = \log (a)$")
fig.tight_layout()

# Plotting derivative and second derivative of g_tilde of x
fig1, ax1 = plt.subplots(2, 1 , figsize=[1.5 * 7.1014, 1.5 * 7.1014 / 1.618])
ax1[0].plot(x, dg_tildedx, label=r"$\tilde{g}'(x)$")
ax1[0].set_ylabel(r"$\frac{d\tilde{g}}{dx}(x)$")
ax1[0].set_xlabel(r"$x = \log (a)$")
ax1[0].set_xlim(-12, 0)
ax1[0].legend()

ax1[1].plot(x, ddg_tildeddx, label=r"$\tilde{g}''(x)$")
ax1[1].legend()
ax1[1].set_xlim(-12, 0)
ax1[1].set_ylabel(r"$\frac{d^2\tilde{g}}{dx^2}(x)$")
ax1[1].set_xlabel(r"$x = \log (a)$")

fig1.tight_layout()
plt.show()


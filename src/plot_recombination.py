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

# Density parameters needed to plot background 
# color corresponding to domination era
cosmo_data = np.loadtxt("cosmology.txt")

OmegaB      = cosmo_data[:, 3]
OmegaCDM    = cosmo_data[:, 4]
OmegaLambda = cosmo_data[:, 5]
OmegaR      = cosmo_data[:, 6]
Omega_sum   = OmegaB + OmegaCDM + OmegaLambda + OmegaR
Omega_m     = OmegaCDM + OmegaB

# Loading data from file
recombo_data = np.loadtxt("recombination.txt")
x               = recombo_data[:, 0]
Xe              = recombo_data[:, 1]
ne              = recombo_data[:, 2] * u.m ** (- 3)
tau             = recombo_data[:, 3]
dtaudx          = recombo_data[:, 4]
ddtaudxdx       = recombo_data[:, 5]
g_tilde         = recombo_data[:, 6]
dg_tildedx      = recombo_data[:, 7]
ddg_tildeddx    = recombo_data[:, 8]
Saha_Xe         = recombo_data[:, 9]


# Computing printout data
g_integral      = np.trapz(g_tilde, x = x)
x_rec           = x[np.argmin(np.abs(Xe - 0.5))]
a_rec           = np.exp(x_rec)
z_rec           = 1 / a_rec - 1
Saha_x_rec      = x[np.argmin(np.abs(Saha_Xe - 0.5))]
Saha_a_rec           = np.exp(Saha_x_rec)
Saha_z_rec           = 1 / Saha_a_rec - 1
log_rel_error   = np.log10(np.abs(g_integral - 1))
x_lss           = x[np.where(g_tilde == g_tilde.max())][0]
a_lss           = np.exp(x_lss)
z_lss           = 1 / a_lss - 1
g_max           = g_tilde.max()
tau_transparent = tau[np.abs(tau - 1).argmin()]
x_transparent   = x[np.abs(tau - 1).argmin()]

# Printout of interesting information
print("----------------------------Some interesting quantities------------------------")
print("Integral of g_tilde: {0}, log rel error: {1}".format(g_integral, log_rel_error))
print("Maximum of g_tilde: x = {0}, g_tilde = {1}".format(x_lss, g_max))
print("Optical depth: tau = {0} at x = {1}".format(tau_transparent, x_transparent))
print("Redshift at last scattering: z = {0}".format(z_lss))
print("Recombination (Xe = 0.5): x = {0}, z = {1}".format(x_rec, z_rec))
print("Recombination Saha (Xe = 0.5): x = {0}, z = {1}".format(Saha_x_rec, Saha_z_rec))
print("-------------------------------------------------------------------------------")

# Generating plots
fig = plt.figure(figsize=[1.5 * 7.1014, 1.5 * 7.1014 / 1.618])

# Plotting electron fraction
ax10 = plt.subplot(221)
ax10.scatter(x_rec, Xe[np.where(x == x_rec)], color = "r", label=r"$(x_{rec}, X_{e,rec})$", zorder = 3)
ax10.plot(x, Xe, label=r"$X_e(x)$", zorder = 1)
ax10.text(-10, Xe[np.where(x == x_rec)], r"$({0:.2f}, {1:.2f})$".format(x_rec, 0.5), color = "r")
ax10.scatter(Saha_x_rec, Saha_Xe[np.where(x == Saha_x_rec)], color = "g", zorder = 4, label=r"$(x_{rec}, X_{e,rec}^{Saha})$")
ax10.plot(x, Saha_Xe, label=r"$X_e^{Saha}(x)$", zorder = 2)
ax10.text(-10, 0.5 * Saha_Xe[np.where(x == Saha_x_rec)], r"$({0:.2f}, {1:.2f})$".format(Saha_x_rec, 0.5), color = "g")
ax10.legend(fontsize = 16)
ax10.set_xlabel(r"$x = \log (a)$")
ax10.set_ylabel(r"$X_e \approx n_e / n_H$")
ax10.set_yscale("log")
ax10.set_xlim(-12, 0)
ax10.set_ylim(1.5e-4, 5)
ax10.axvspan(
    np.min(x),
    x[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    alpha=0.3,
    color="orange",
)
ax10.axvspan(
    x[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    x[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    alpha=0.3,
    color="cyan",
)
ax10.axvspan(
    x[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    np.max(x),
    alpha=0.3,
    color="purple",
)


# Plotting electron density
ax11 = plt.subplot(222)
ax11.plot(x, ne, label=r"$n_e(x)$")
ax11.legend()
ax11.set_xlabel(r"$x = \log (a)$")
ax11.set_ylabel(r"$n_e [\mathrm{m^{-3}}]$")
ax11.set_yscale("log")
ax11.set_xlim(-12, 0)
ax11.axvspan(
    np.min(x),
    x[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    alpha=0.3,
    color="orange",
)
ax11.axvspan(
    x[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    x[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    alpha=0.3,
    color="cyan",
)
ax11.axvspan(
    x[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    np.max(x),
    alpha=0.3,
    color="purple",
)


# Plotting optical depth
ax12 = plt.subplot(212)
ax12.plot(x, tau, label=r"$\tau(x)$")
ax12.plot(x, - dtaudx, label=r"$-\tau'(x)$", linestyle = "--")
ax12.plot(x, ddtaudxdx, label=r"$\tau''(x)$", linestyle = "-.")
ax12.scatter(x_transparent, tau_transparent, color = "r", label = r"$(x_{lss}, \tau_{lss})$")
ax12.text(-8.2, 1, "({0:.2f}, 1)".format(x_transparent), color = "r")
ax12.legend()
ax12.set_xlabel(r"$x = \log (a)$")
ax12.set_ylabel(r"$\tau(x)$")
ax12.set_ylim(1e-8, 1e8)
ax12.set_xlim(-12, 0)
ax12.set_yscale("log")
ax12.axvspan(
    np.min(x),
    x[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    alpha=0.3,
    color="orange",
)
ax12.axvspan(
    x[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    x[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    alpha=0.3,
    color="cyan",
)
ax12.axvspan(
    x[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    np.max(x),
    alpha=0.3,
    color="purple",
)

fig.tight_layout()
fig.savefig("../doc/Figures/Xe_ne_tau.pdf", dpi=1000)


# Plotting visibility function, derivative and second derivative thereof
fig1, ax1 = plt.subplots(2, 2 , figsize=[1.5 * 7.1014, 1.5 * 7.1014 / 1.618])

ax1[0, 0].scatter(x_lss, g_max, color = "r", label = r"$(x_{lss}, \tilde{g}_{lss})$")
ax1[0, 0].plot(x, g_tilde, label=r"$\tilde{g}(x)$")
ax1[0, 0].text(-10, g_max, r"$({0:.2f}, {1:.2f})$".format(x_lss, g_max), color = "r")
ax1[0, 0].set_xlim(-12, 0)
ax1[0, 0].legend()
ax1[0, 0].set_ylabel(r"$\tilde{g}(x)$")
ax1[0, 0].set_xlabel(r"$x = \log (a)$")
ax1[0, 0].axvspan(
    np.min(x),
    x[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    alpha=0.3,
    color="orange",
)
ax1[0, 0].axvspan(
    x[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    x[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    alpha=0.3,
    color="cyan",
)
ax1[0, 0].axvspan(
    x[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    np.max(x),
    alpha=0.3,
    color="purple",
)

ax1[0, 1].plot(x, dg_tildedx, label=r"$\tilde{g}'(x)$")
ax1[0, 1].set_ylabel(r"$\frac{d\tilde{g}}{dx}(x)$")
ax1[0, 1].set_xlabel(r"$x = \log (a)$")
ax1[0, 1].set_xlim(-12, 0)
ax1[0, 1].legend()
ax1[0, 1].axvspan(
    np.min(x),
    x[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    alpha=0.3,
    color="orange",
)
ax1[0, 1].axvspan(
    x[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    x[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    alpha=0.3,
    color="cyan",
)
ax1[0, 1].axvspan(
    x[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    np.max(x),
    alpha=0.3,
    color="purple",
)


ax1[1, 0].plot(x, ddg_tildeddx, label=r"$\tilde{g}''(x)$")
ax1[1, 0].legend()
ax1[1, 0].set_xlim(-12, 0)
ax1[1, 0].set_ylabel(r"$\frac{d^2\tilde{g}}{dx^2}(x)$")
ax1[1, 0].set_xlabel(r"$x = \log (a)$")
ax1[1, 0].axvspan(
    np.min(x),
    x[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    alpha=0.3,
    color="orange",
)
ax1[1, 0].axvspan(
    x[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    x[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    alpha=0.3,
    color="cyan",
)
ax1[1, 0].axvspan(
    x[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    np.max(x),
    alpha=0.3,
    color="purple",
)

# Peak normalized visibility function and its derivatives
ax1[1, 1].plot(x, g_tilde / np.max(g_max), label=r"$\tilde{g}(x)$")
ax1[1, 1].plot(x, dg_tildedx / np.max(np.max(dg_tildedx)), label=r"$\tilde{g}'(x)$")
ax1[1, 1].plot(x, ddg_tildeddx / np.max(np.abs(ddg_tildeddx)), label=r"$\tilde{g}''(x)$")
ax1[1, 1].legend(loc = 1)
ax1[1, 1].set_xlim(-7.5, -5.7)
ax1[1, 1].set_ylabel(r"Peak normalized")
ax1[1, 1].set_xlabel(r"$x = \log (a)$")
ax1[1, 1].axvspan(
    np.min(x),
    x[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    alpha=0.3,
    color="orange",
)
ax1[1, 1].axvspan(
    x[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    x[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    alpha=0.3,
    color="cyan",
)
ax1[1, 1].axvspan(
    x[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    np.max(x),
    alpha=0.3,
    color="purple",
)

fig1.tight_layout()
fig1.savefig("../doc/Figures/g_tilde.pdf", dpi=1000)
plt.show()

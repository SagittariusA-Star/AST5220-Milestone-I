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

# Loading data from file
pertub_data = np.loadtxt("perturbations_k0.01.txt")
k      = 0.01 / u.Mpc
x      = pertub_data[:, 0]
Theta0 = pertub_data[:, 1]
Theta1 = pertub_data[:, 2]
Theta2 = pertub_data[:, 3]
Phi    = pertub_data[:, 4]
Psi    = pertub_data[:, 5]
delta_cdm = pertub_data[:, 6]
delta_b   = pertub_data[:, 7]
v_cdm     = pertub_data[:, 8] 
v_b       = pertub_data[:, 9] 

# Plotting visibility function, derivative and second derivative thereof
fig, ax = plt.subplots(2, 2 , figsize=[1.5 * 7.1014, 1.5 * 7.1014 / 1.618])

ax[0, 0].plot(x, Theta0, label=r"$\Theta_0$")
ax[0, 0].set_xlim(-12, 0)
ax[0, 0].legend()
ax[0, 0].set_ylabel(r"$\Theta_0$")
ax[0, 0].set_xlabel(r"$x = \log (a)$")
"""
ax[0, 0].axvspan(
    np.min(x),
    x[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    alpha=0.3,
    color="orange",
)
ax[0, 0].axvspan(
    x[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    x[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    alpha=0.3,
    color="cyan",
)
ax[0, 0].axvspan(
    x[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    np.max(x),
    alpha=0.3,
    color="purple",
)
"""
ax[0, 1].plot(x, Theta1, label=r"$\Theta_1$")
ax[0, 1].set_ylabel(r"$\Theta_1$")
ax[0, 1].set_xlabel(r"$x = \log (a)$")
ax[0, 1].set_xlim(-12, 0)
ax[0, 1].legend()
"""
ax[0, 1].axvspan(
    np.min(x),
    x[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    alpha=0.3,
    color="orange",
)
ax[0, 1].axvspan(
    x[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    x[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    alpha=0.3,
    color="cyan",
)
ax[0, 1].axvspan(
    x[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    np.max(x),
    alpha=0.3,
    color="purple",
)
"""

ax[1, 0].plot(x, Theta2, label=r"$\Theta_2$")
ax[1, 0].legend()
ax[1, 0].set_xlim(-12, 0)
ax[1, 0].set_ylabel(r"$\Theta_2$")
ax[1, 0].set_xlabel(r"$x = \log (a)$")
"""
ax[1, 0].axvspan(
    np.min(x),
    x[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    alpha=0.3,
    color="orange",
)
ax[1, 0].axvspan(
    x[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    x[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    alpha=0.3,
    color="cyan",
)
ax[1, 0].axvspan(
    x[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    np.max(x),
    alpha=0.3,
    color="purple",
)
"""
# Peak normalized visibility function and its derivatives
ax[1, 1].plot(x, Phi, label=r"$\Phi$")
#ax[1, 1].plot(x, Psi, label=r"$\Psi$")
ax[1, 1].legend()
ax[1, 1].set_xlim(-12, 0)
ax[1, 1].set_ylabel(r"$\Phi$")
ax[1, 1].set_xlabel(r"$x = \log (a)$")
"""
ax[1, 1].axvspan(
    np.min(x),
    x[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    alpha=0.3,
    color="orange",
)
ax[1, 1].axvspan(
    x[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    x[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    alpha=0.3,
    color="cyan",
)
ax[1, 1].axvspan(
    x[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    np.max(x),
    alpha=0.3,
    color="purple",
)
"""
fig.tight_layout()

plt.figure()
plt.semilogy(x, delta_cdm, label = r"$\delta_{cdm}$")
plt.semilogy(x, delta_b, label = r"$\delta_{b}$", linestyle = "--")
plt.legend(loc = 0)

plt.figure()
plt.semilogy(x, v_cdm, label = r"$v_{cdm}$")
plt.semilogy(x, v_b, label = r"$v_{b}$", linestyle = "--")
plt.legend(loc = 0)

plt.show()

"""
# Generating plots
fig = plt.figure(figsize=[1.5 * 7.1014, 1.5 * 7.1014 / 1.618])

# Plotting electron fraction
ax0 = plt.subplot(221)
ax0.scatter(x_rec, Xe[np.where(x == x_rec)], color = "r", label=r"$(x_{rec}, X_{e,rec})$", zorder = 3)
ax0.plot(x, Xe, label=r"$X_e(x)$", zorder = 1)
ax0.text(-10, Xe[np.where(x == x_rec)], r"$({0:.2f}, {1:.2f})$".format(x_rec, 0.5), color = "r")
ax0.scatter(Saha_x_rec, Saha_Xe[np.where(x == Saha_x_rec)], color = "g", zorder = 4, label=r"$(x_{rec}, X_{e,rec}^{Saha})$")
ax0.plot(x, Saha_Xe, label=r"$X_e^{Saha}(x)$", zorder = 2)
ax0.text(-10, 0.5 * Saha_Xe[np.where(x == Saha_x_rec)], r"$({0:.2f}, {1:.2f})$".format(Saha_x_rec, 0.5), color = "g")
ax0.legend(fontsize = 16)
ax0.set_xlabel(r"$x = \log (a)$")
ax0.set_ylabel(r"$X_e \approx n_e / n_H$")
ax0.set_yscale("log")
ax0.set_xlim(-12, 0)
ax0.set_ylim(1.5e-4, 5)
ax0.axvspan(
    np.min(x),
    x[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    alpha=0.3,
    color="orange",
)
ax0.axvspan(
    x[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    x[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    alpha=0.3,
    color="cyan",
)
ax0.axvspan(
    x[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    np.max(x),
    alpha=0.3,
    color="purple",
)


# Plotting electron density
ax1 = plt.subplot(222)
ax1.plot(x, ne, label=r"$n_e(x)$")
ax1.legend()
ax1.set_xlabel(r"$x = \log (a)$")
ax1.set_ylabel(r"$n_e [\mathrm{m^{-3}}]$")
ax1.set_yscale("log")
ax1.set_xlim(-12, 0)
ax1.axvspan(
    np.min(x),
    x[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    alpha=0.3,
    color="orange",
)
ax1.axvspan(
    x[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    x[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    alpha=0.3,
    color="cyan",
)
ax1.axvspan(
    x[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    np.max(x),
    alpha=0.3,
    color="purple",
)


# Plotting optical depth
ax2 = plt.subplot(212)
ax2.plot(x, tau, label=r"$\tau(x)$")
ax2.plot(x, - dtaudx, label=r"$-\tau'(x)$", linestyle = "--")
ax2.plot(x, ddtaudxdx, label=r"$\tau''(x)$", linestyle = "-.")
ax2.scatter(x_transparent, tau_transparent, color = "r", label = r"$(x_{lss}, \tau_{lss})$")
ax2.text(-8.2, 1, "({0:.2f}, 1)".format(x_transparent), color = "r")
ax2.legend()
ax2.set_xlabel(r"$x = \log (a)$")
ax2.set_ylabel(r"$\tau(x)$")
ax2.set_ylim(1e-8, 1e8)
ax2.set_xlim(-12, 0)
ax2.set_yscale("log")
ax2.axvspan(
    np.min(x),
    x[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    alpha=0.3,
    color="orange",
)
ax2.axvspan(
    x[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    x[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    alpha=0.3,
    color="cyan",
)
ax2.axvspan(
    x[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    np.max(x),
    alpha=0.3,
    color="purple",
)

fig.tight_layout()
fig.savefig("../doc/Figures/Xe_ne_tau.pdf", dpi=1000)


"""
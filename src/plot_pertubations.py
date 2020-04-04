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

# Density parameters needed to plot background 
# color corresponding to domination era
cosmo_data = np.loadtxt("cosmology.txt")

x_cosmo      = cosmo_data[:, 0]
OmegaB      = cosmo_data[:, 3]
OmegaCDM    = cosmo_data[:, 4]
OmegaLambda = cosmo_data[:, 5]
OmegaR      = cosmo_data[:, 6]
Omega_sum   = OmegaB + OmegaCDM + OmegaLambda + OmegaR
Omega_m     = OmegaCDM + OmegaB

# Loading data from file
recombo_data = np.loadtxt("recombination.txt")
x_recombo       = recombo_data[:, 0]
Xe              = recombo_data[:, 1]
ne              = recombo_data[:, 2] * u.m ** (- 3)
tau             = recombo_data[:, 3]
dtaudx          = recombo_data[:, 4]
ddtaudxdx       = recombo_data[:, 5]
g_tilde         = recombo_data[:, 6]
dg_tildedx      = recombo_data[:, 7]
ddg_tildeddx    = recombo_data[:, 8]
Saha_Xe         = recombo_data[:, 9]

recombo_index = np.where(g_tilde >= 1.5e-1)
Xe_half_index = np.argmin(np.abs(Xe - 0.5))

# Loading data from file
pertub_data_001 = np.loadtxt("perturbations_k0.01.txt")
pertub_data_0001 = np.loadtxt("perturbations_k0.001.txt")
pertub_data_01 = np.loadtxt("perturbations_k0.1.txt")

k_001      = 0.01 / u.Mpc
x_001      = pertub_data_001[:, 0]
Theta0_001 = pertub_data_001[:, 1]
Theta1_001 = pertub_data_001[:, 2]
Theta2_001 = pertub_data_001[:, 3]
Phi_001    = pertub_data_001[:, 4]
Psi_001    = pertub_data_001[:, 5]
delta_cdm_001 = pertub_data_001[:, 6]
delta_b_001   = pertub_data_001[:, 7]
v_cdm_001     = pertub_data_001[:, 8] 
v_b_001       = pertub_data_001[:, 9] 

k_0001      = 0.001 / u.Mpc
x_0001      = pertub_data_0001[:, 0]
Theta0_0001 = pertub_data_0001[:, 1]
Theta1_0001 = pertub_data_0001[:, 2]
Theta2_0001 = pertub_data_0001[:, 3]
Phi_0001    = pertub_data_0001[:, 4]
Psi_0001    = pertub_data_0001[:, 5]
delta_cdm_0001 = pertub_data_0001[:, 6]
delta_b_0001   = pertub_data_0001[:, 7]
v_cdm_0001     = pertub_data_0001[:, 8] 
v_b_0001       = pertub_data_0001[:, 9] 

k_01      = 0.1 / u.Mpc
x_01      = pertub_data_01[:, 0]
Theta0_01 = pertub_data_01[:, 1]
Theta1_01 = pertub_data_01[:, 2]
Theta2_01 = pertub_data_01[:, 3]
Phi_01    = pertub_data_01[:, 4]
Psi_01    = pertub_data_01[:, 5]
delta_cdm_01 = pertub_data_01[:, 6]
delta_b_01   = pertub_data_01[:, 7]
v_cdm_01     = pertub_data_01[:, 8] 
v_b_01       = pertub_data_01[:, 9] 


# Plotting visibility function, derivative and second derivative thereof
fig, ax = plt.subplots(2, 2 , figsize=[1.5 * 7.1014, 1.5 * 7.1014 / 1.618])
fig.suptitle(r"$k = 0.01 \rm{Mpc}^{-1}$")

ax[0, 0].plot(x_001, Theta0_001, label=rf"$\Theta_0(k = {k_001.value}/\mathrm{{Mpc}})$")
ax[0, 0].plot(x_0001, Theta0_0001, label=rf"$\Theta_0(k = {k_0001.value}/\mathrm{{Mpc}})$")
ax[0, 0].plot(x_01, Theta0_01, label=rf"$\Theta_0(k = {k_01.value}/\mathrm{{Mpc}})$")
ax[0, 0].set_xlim(-12, 0)
ax[0, 0].legend()
ax[0, 0].set_ylabel(r"$\Theta_0$")
ax[0, 0].set_xlabel(r"$x = \log (a)$")
"""
ax[0, 0].axvspan(
    np.min(x_cosmo),
    x_cosmo[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    alpha=0.3,
    color="orange",
)
ax[0, 0].axvspan(
    x_cosmo[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    x_cosmo[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    alpha=0.3,
    color="cyan",
)
ax[0, 0].axvspan(
    x_cosmo[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    np.max(x_cosmo),
    alpha=0.3,
    color="purple",
)

ax[0, 0].axvspan(
    x_cosmo[recombo_index][0],
    x_cosmo[recombo_index][-1],
    alpha=0.4,
    color="red",
)

ax[0, 0].axvline(x = x_cosmo[Xe_half_index], ymin = -1e5, ymax = 1e5, color = "orangered", linestyle = ":")
"""
ax[0, 1].plot(x_001, Theta1_001, label=rf"$\Theta_1(k = {k_001.value}/\mathrm{{Mpc}})$")
ax[0, 1].plot(x_0001, Theta1_0001, label=rf"$\Theta_1(k = {k_0001.value}/\mathrm{{Mpc}})$")
ax[0, 1].plot(x_01, Theta1_01, label=rf"$\Theta_1(k = {k_01.value}/\mathrm{{Mpc}})$")
ax[0, 1].set_ylabel(r"$\Theta_1$")
ax[0, 1].set_xlabel(r"$x = \log (a)$")
ax[0, 1].set_xlim(-12, 0)
ax[0, 1].legend()
"""
ax[0, 1].axvspan(
    np.min(x_cosmo),
    x_cosmo[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    alpha=0.3,
    color="orange",
)
ax[0, 1].axvspan(
    x_cosmo[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    x_cosmo[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    alpha=0.3,
    color="cyan",
)
ax[0, 1].axvspan(
    x_cosmo[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    np.max(x_cosmo),
    alpha=0.3,
    color="purple",
)
ax[0, 1].axvspan(
    x_cosmo[recombo_index][0],
    x_cosmo[recombo_index][-1],
    alpha=0.4,
    color="red",
)

ax[0, 1].axvline(x = x_cosmo[Xe_half_index], ymin = -1e5, ymax = 1e5, color = "orangered", linestyle = ":")
"""
ax[1, 0].plot(x_001, Theta2_001, label=rf"$\Theta_2(k = {k_001.value}/\mathrm{{Mpc}})$")
ax[1, 0].plot(x_0001, Theta2_0001, label=rf"$\Theta_2(k = {k_0001.value}/\mathrm{{Mpc}})$")
ax[1, 0].plot(x_01, Theta2_01, label=rf"$\Theta_2(k = {k_01.value}/\mathrm{{Mpc}})$")
ax[1, 0].legend()
ax[1, 0].set_xlim(-12, 0)
ax[1, 0].set_ylabel(r"$\Theta_2$")
ax[1, 0].set_xlabel(r"$x = \log (a)$")
"""
ax[1, 0].axvspan(
    np.min(x_cosmo),
    x_cosmo[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    alpha=0.3,
    color="orange",
)
ax[1, 0].axvspan(
    x_cosmo[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    x_cosmo[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    alpha=0.3,
    color="cyan",
)
ax[1, 0].axvspan(
    x_cosmo[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    np.max(x_cosmo),
    alpha=0.3,
    color="purple",
)
ax[1, 0].axvspan(
    x_cosmo[recombo_index][0],
    x_cosmo[recombo_index][-1],
    alpha=0.4,
    color="red",
)

ax[1, 0].axvline(x = x_cosmo[Xe_half_index], ymin = -1e5, ymax = 1e5, color = "orangered", linestyle = ":")
"""
# Peak normalized visibility function and its derivatives
ax[1, 1].plot(x_001, Phi_001, label=rf"$\Phi(k = {k_001.value}/\mathrm{{Mpc}})$")
ax[1, 1].plot(x_0001, Phi_0001, label=rf"$\Phi(k = {k_0001.value}/\mathrm{{Mpc}})$")
ax[1, 1].plot(x_01, Phi_01, label=rf"$\Phi(k = {k_01.value}/\mathrm{{Mpc}})$")

#ax[1, 1].plot(x_001, Psi_001, label=rf"$\Psi(k = {k_001.value}/\mathrm{{Mpc}})$")
#ax[1, 1].plot(x_0001, Psi_0001, label=rf"$\Psi(k = {k_0001.value}/\mathrm{{Mpc}})$")
#ax[1, 1].plot(x_01, Psi_01, label=rf"$\Psi(k = {k_01.value}/\mathrm{{Mpc}})$")
ax[1, 1].legend()
#ax[1, 1].set_xlim(-12, 0)
ax[1, 1].set_ylabel(r"$\Phi$")
ax[1, 1].set_xlabel(r"$x = \log (a)$")
"""
ax[1, 1].axvspan(
    np.min(x_cosmo),
    x_cosmo[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    alpha=0.3,
    color="orange",
)
ax[1, 1].axvspan(
    x_cosmo[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    x_cosmo[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    alpha=0.3,
    color="cyan",
)
ax[1, 1].axvspan(
    x_cosmo[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    np.max(x_cosmo),
    alpha=0.3,
    color="purple",
)

ax[1, 1].axvspan(
    x_cosmo[recombo_index][0],
    x_cosmo[recombo_index][-1],
    alpha=0.4,
    color="red",
)

ax[1, 1].axvline(x = x_cosmo[Xe_half_index], ymin = -1e5, ymax = 1e5, color = "orangered", linestyle = ":")
"""
fig.tight_layout()
fig.savefig("../doc/Figures/fig1.pdf", dpi=1000)


fig1, ax1 = plt.subplots(1, 2 , figsize=[1.5 * 7.1014, 1.5 * 7.1014 / 1.618])
fig1.suptitle(r"$k = 0.01 \rm{\mathrm{{Mpc}}}^{-1}$")

ax1[0].plot(x_001, np.abs(delta_cdm_001), label = fr"$\delta_{{cdm}}(k = {k_001.value}/\mathrm{{Mpc}})$")
ax1[0].plot(x_0001, np.abs(delta_cdm_0001), label = fr"$\delta_{{cdm}}(k = {k_0001.value}/\mathrm{{Mpc}})$")
ax1[0].plot(x_01, np.abs(delta_cdm_01), label = fr"$\delta_{{cdm}}(k = {k_01.value}/\mathrm{{Mpc}})$")

ax1[0].plot(x_001, np.abs(delta_b_001), label = rf"$\delta_{{b}}(k = {k_001.value}/\mathrm{{Mpc}})$", linestyle = "--")
ax1[0].plot(x_0001, np.abs(delta_b_0001), label = rf"$\delta_{{b}}(k = {k_0001.value}/\mathrm{{Mpc}})$", linestyle = "--")
ax1[0].plot(x_01, np.abs(delta_b_01), label = rf"$\delta_{{b}}(k = {k_01.value}/\mathrm{{Mpc}})$", linestyle = "--")
#ax1[0].set_xlim(-12, 0)
ax1[0].legend()
ax1[0].set_ylabel(r"$\delta$")
ax1[0].set_xlabel(r"$x = \log (a)$")
ax1[0].set_yscale("log")
"""
ax1[0].axvspan(
    np.min(x_cosmo),
    x_cosmo[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    alpha=0.3,
    color="orange",
)
ax1[0].axvspan(
    x_cosmo[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    x_cosmo[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    alpha=0.3,
    color="cyan",
)
ax1[0].axvspan(
    x_cosmo[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    np.max(x_cosmo),
    alpha=0.3,
    color="purple",
)

ax1[0].axvspan(
    x_cosmo[recombo_index][0],
    x_cosmo[recombo_index][-1],
    alpha=0.4,
    color="red",
)
ax1[0].axvline(x = x_cosmo[Xe_half_index], ymin = 1e-5, ymax = 1e5, color = "orangered", linestyle = ":")
"""
ax1[1].plot(x_001, np.abs(v_cdm_001), label = rf"$v_{{cdm}}(k = {k_001.value}/\mathrm{{Mpc}})$")
ax1[1].plot(x_0001, np.abs(v_cdm_0001), label = rf"$v_{{cdm}}(k = {k_0001.value}/\mathrm{{Mpc}})$")
ax1[1].plot(x_01, np.abs(v_cdm_01), label = rf"$v_{{cdm}}(k = {k_01.value}/\mathrm{{Mpc}})$")

ax1[1].plot(x_001, np.abs(v_b_001), label = rf"$v_{{b}}(k = {k_001.value}/\mathrm{{Mpc}})$", linestyle = "--")
ax1[1].plot(x_0001, np.abs(v_b_0001), label = rf"$v_{{b}}(k = {k_0001.value}/\mathrm{{Mpc}})$", linestyle = "--")
ax1[1].plot(x_01, np.abs(v_b_01), label = rf"$v_{{b}}(k = {k_01.value}/\mathrm{{Mpc}})$", linestyle = "--")
#ax1[1].set_xlim(-12, 0)
ax1[1].legend()
ax1[1].set_ylabel(r"$v$")
ax1[1].set_xlabel(r"$x = \log (a)$")
ax1[1].set_yscale("log")
fig1.savefig("../doc/Figures/fig2.pdf", dpi=1000)
"""
ax1[1].axvspan(
    np.min(x_cosmo),
    x_cosmo[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    alpha=0.3,
    color="orange",
)
ax1[1].axvspan(
    x_cosmo[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    x_cosmo[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    alpha=0.3,
    color="cyan",
)
ax1[1].axvspan(
    x_cosmo[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    np.max(x_cosmo),
    alpha=0.3,
    color="purple",
)
ax1[1].axvspan(
    x_cosmo[recombo_index][0],
    x_cosmo[recombo_index][-1],
    alpha=0.4,
    color="red",
)
ax1[1].axvline(x = x_cosmo[Xe_half_index], ymin = 1e-5, ymax = 1e5, color = "orangered", linestyle = ":")
"""
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
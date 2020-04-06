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
#pertub_data_001 = np.loadtxt("perturbations_k0.01.txt")
#pertub_data_0001 = np.loadtxt("perturbations_k0.001.txt")
#pertub_data_01 = np.loadtxt("perturbations_k0.1.txt")
pertub_data_0 = np.loadtxt("perturbations_k5e-5.txt")
pertub_data_1 = np.loadtxt("perturbations_k2.85e-4.txt")
pertub_data_2 = np.loadtxt("perturbations_k1.623e-3.txt")
pertub_data_3 = np.loadtxt("perturbations_k9.244e-3.txt")
pertub_data_4 = np.loadtxt("perturbations_k5.266e-2.txt")
pertub_data_5 = np.loadtxt("perturbations_k3e-1.txt")

"""
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
"""
k_0      = 5e-5 / u.Mpc
x_0      = pertub_data_0[:, 0]
Theta0_0 = pertub_data_0[:, 1]
Theta1_0 = pertub_data_0[:, 2]
Theta2_0 = pertub_data_0[:, 3]
Phi_0    = pertub_data_0[:, 4]
Psi_0    = pertub_data_0[:, 5]
delta_cdm_0 = pertub_data_0[:, 6]
delta_b_0   = pertub_data_0[:, 7]
v_cdm_0     = pertub_data_0[:, 8] 
v_b_0       = pertub_data_0[:, 9] 

k_1      = 2.85e-4 / u.Mpc
x_1      = pertub_data_1[:, 0]
Theta0_1 = pertub_data_1[:, 1]
Theta1_1 = pertub_data_1[:, 2]
Theta2_1 = pertub_data_1[:, 3]
Phi_1    = pertub_data_1[:, 4]
Psi_1    = pertub_data_1[:, 5]
delta_cdm_1 = pertub_data_1[:, 6]
delta_b_1   = pertub_data_1[:, 7]
v_cdm_1     = pertub_data_1[:, 8] 
v_b_1       = pertub_data_1[:, 9] 

k_2      = 1.623e-3 / u.Mpc
x_2      = pertub_data_2[:, 0]
Theta0_2 = pertub_data_2[:, 1]
Theta1_2 = pertub_data_2[:, 2]
Theta2_2 = pertub_data_2[:, 3]
Phi_2    = pertub_data_2[:, 4]
Psi_2    = pertub_data_2[:, 5]
delta_cdm_2 = pertub_data_2[:, 6]
delta_b_2   = pertub_data_2[:, 7]
v_cdm_2     = pertub_data_2[:, 8] 
v_b_2       = pertub_data_2[:, 9] 

k_3      = 9.244e-3 / u.Mpc
x_3      = pertub_data_3[:, 0]
Theta0_3 = pertub_data_3[:, 1]
Theta1_3 = pertub_data_3[:, 2]
Theta2_3 = pertub_data_3[:, 3]
Phi_3    = pertub_data_3[:, 4]
Psi_3    = pertub_data_3[:, 5]
delta_cdm_3 = pertub_data_3[:, 6]
delta_b_3   = pertub_data_3[:, 7]
v_cdm_3     = pertub_data_3[:, 8] 
v_b_3       = pertub_data_3[:, 9] 

k_4      = 5.266e-2 / u.Mpc
x_4      = pertub_data_4[:, 0]
Theta0_4 = pertub_data_4[:, 1]
Theta1_4 = pertub_data_4[:, 2]
Theta2_4 = pertub_data_4[:, 3]
Phi_4    = pertub_data_4[:, 4]
Psi_4    = pertub_data_4[:, 5]
delta_cdm_4 = pertub_data_4[:, 6]
delta_b_4   = pertub_data_4[:, 7]
v_cdm_4     = pertub_data_4[:, 8] 
v_b_4       = pertub_data_4[:, 9] 

k_5      = 3e-1 / u.Mpc
x_5      = pertub_data_5[:, 0]
Theta0_5 = pertub_data_5[:, 1]
Theta1_5 = pertub_data_5[:, 2]
Theta2_5 = pertub_data_5[:, 3]
Phi_5    = pertub_data_5[:, 4]
Psi_5    = pertub_data_5[:, 5]
delta_cdm_5 = pertub_data_5[:, 6]
delta_b_5   = pertub_data_5[:, 7]
v_cdm_5     = pertub_data_5[:, 8] 
v_b_5       = pertub_data_5[:, 9] 

# Plotting visibility function, derivative and second derivative thereof
fig, ax = plt.subplots(2, 2 , figsize=[1.5 * 7.1014, 1.5 * 7.1014 / 1.618])

ax[0, 0].plot(x_5, Theta0_5, label=rf"$k = {k_5.value:.2e}/\mathrm{{Mpc}}$", color = "orange")
ax[0, 0].plot(x_4, Theta0_4, label=rf"$k = {k_4.value:.2e}/\mathrm{{Mpc}}$", color = "k")
ax[0, 0].plot(x_3, Theta0_3, label=rf"$k = {k_3.value:.2e}/\mathrm{{Mpc}}$", color = "m")
ax[0, 0].plot(x_2, Theta0_2, label=rf"$k = {k_2.value:.2e}/\mathrm{{Mpc}}$", color = "b")
ax[0, 0].plot(x_1, Theta0_1, label=rf"$k = {k_1.value:.2e}/\mathrm{{Mpc}}$", color = "r")
ax[0, 0].plot(x_0, Theta0_0, label=rf"$k = {k_0.value:.2e}/\mathrm{{Mpc}}$", color = "g")
ax[0, 0].set_xlim(-12, 0)
#ax[0, 0].legend()
ax[0, 0].set_ylabel(r"$\Theta_0$")
ax[0, 0].set_xlabel(r"$x = \log (a)$")
fig.legend(loc='center', bbox_to_anchor=(0.3, 0.3), fontsize = 14)
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
ax[0, 1].plot(x_5, Theta1_5, label=rf"$k = {k_5.value:.2e}/\mathrm{{Mpc}}$", color = "orange")
ax[0, 1].plot(x_4, Theta1_4, label=rf"$k = {k_4.value:.2e}/\mathrm{{Mpc}}$", color = "k")
ax[0, 1].plot(x_3, Theta1_3, label=rf"$k = {k_3.value:.2e}/\mathrm{{Mpc}}$", color = "m")
ax[0, 1].plot(x_2, Theta1_2, label=rf"$k = {k_2.value:.2e}/\mathrm{{Mpc}}$", color = "b")
ax[0, 1].plot(x_1, Theta1_1, label=rf"$k = {k_1.value:.2e}/\mathrm{{Mpc}}$", color = "r")
ax[0, 1].plot(x_0, Theta1_0, label=rf"$k = {k_0.value:.2e}/\mathrm{{Mpc}}$", color = "g")
ax[0, 1].set_ylabel(r"$\Theta_1$")
ax[0, 1].set_xlabel(r"$x = \log (a)$")
ax[0, 1].set_xlim(-12, 0)
#ax[0, 1].legend()
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
ax[1, 0].set_frame_on(False)
ax[1, 0].get_xaxis().set_visible(False)
ax[1, 0].get_yaxis().set_visible(False)
"""
ax[1, 0].plot(x_5, Theta2_5, label=rf"$(k = {k_5.value:.2e}/\mathrm{{Mpc}}$", color = "orange")
ax[1, 0].plot(x_4, Theta2_4, label=rf"$(k = {k_4.value:.2e}/\mathrm{{Mpc}}$", color = "k")
ax[1, 0].plot(x_3, Theta2_3, label=rf"$(k = {k_3.value:.2e}/\mathrm{{Mpc}}$", color = "m")
ax[1, 0].plot(x_2, Theta2_2, label=rf"$(k = {k_2.value:.2e}/\mathrm{{Mpc}}$", color = "b")
ax[1, 0].plot(x_1, Theta2_1, label=rf"$(k = {k_1.value:.2e}/\mathrm{{Mpc}}$", color = "r")
ax[1, 0].plot(x_0, Theta2_0, label=rf"$(k = {k_0.value:.2e}/\mathrm{{Mpc}}$", color = "g")
#ax[1, 0].legend()
ax[1, 0].set_xlim(-12, 0)
ax[1, 0].set_ylabel(r"$\Theta_2$")
ax[1, 0].set_xlabel(r"$x = \log (a)$")
"""
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
ax[1, 1].plot(x_5, Phi_5, label=rf"$k = {k_5.value:.2e}/\mathrm{{Mpc}}$", color = "orange")
ax[1, 1].plot(x_4, Phi_4, label=rf"$k = {k_4.value:.2e}/\mathrm{{Mpc}}$", color = "k")
ax[1, 1].plot(x_3, Phi_3, label=rf"$k = {k_3.value:.2e}/\mathrm{{Mpc}}$", color = "m")
ax[1, 1].plot(x_2, Phi_2, label=rf"$k = {k_2.value:.2e}/\mathrm{{Mpc}}$", color = "b")
ax[1, 1].plot(x_1, Phi_1, label=rf"$k = {k_1.value:.2e}/\mathrm{{Mpc}}$", color = "r")
ax[1, 1].plot(x_0, Phi_0, label=rf"$k = {k_0.value:.2e}/\mathrm{{Mpc}}$", color = "g")

ax[1, 1].plot(x_5, Psi_5, linestyle = "--", color = "orange")
ax[1, 1].plot(x_4, Psi_4, linestyle = "--", color = "k")
ax[1, 1].plot(x_3, Psi_3, linestyle = "--", color = "m")
ax[1, 1].plot(x_2, Psi_2, linestyle = "--", color = "b")
ax[1, 1].plot(x_1, Psi_1, linestyle = "--", color = "r")
ax[1, 1].plot(x_0, Psi_0, linestyle = "--", color = "g")
#ax[1, 1].legend()
#ax[1, 1].set_xlim(-12, 0)
ax[1, 1].set_ylabel(r"$\Phi$  $(\Psi)$")
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

ax1[0].plot(x_5, np.abs(delta_cdm_5), label = fr"$k = {k_5.value:.2e}/\mathrm{{Mpc}}$", color = "orange")
ax1[0].plot(x_4, np.abs(delta_cdm_4), label = fr"$k = {k_4.value:.2e}/\mathrm{{Mpc}}$", color = "k")
ax1[0].plot(x_3, np.abs(delta_cdm_3), label = fr"$k = {k_3.value:.2e}/\mathrm{{Mpc}}$", color = "m")
ax1[0].plot(x_2, np.abs(delta_cdm_2), label = fr"$k = {k_2.value:.2e}/\mathrm{{Mpc}}$", color = "b")
ax1[0].plot(x_1, np.abs(delta_cdm_1), label = fr"$k = {k_1.value:.2e}/\mathrm{{Mpc}}$", color = "r")
ax1[0].plot(x_0, np.abs(delta_cdm_0), label = fr"$k = {k_0.value:.2e}/\mathrm{{Mpc}}$", color = "g")

ax1[0].plot(x_5, np.abs(delta_b_5), linestyle = "--", color = "orange")
ax1[0].plot(x_4, np.abs(delta_b_4), linestyle = "--", color = "k")
ax1[0].plot(x_3, np.abs(delta_b_3), linestyle = "--", color = "m")
ax1[0].plot(x_2, np.abs(delta_b_2), linestyle = "--", color = "b")
ax1[0].plot(x_1, np.abs(delta_b_1), linestyle = "--", color = "r")
ax1[0].plot(x_0, np.abs(delta_b_0), linestyle = "--", color = "g")
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
ax1[1].plot(x_5, np.abs(v_cdm_5), label = rf"$k = {k_5.value:.2e}/\mathrm{{Mpc}}$", color = "orange")
ax1[1].plot(x_4, np.abs(v_cdm_4), label = rf"$k = {k_4.value:.2e}/\mathrm{{Mpc}}$", color = "k")
ax1[1].plot(x_3, np.abs(v_cdm_3), label = rf"$k = {k_3.value:.2e}/\mathrm{{Mpc}}$", color = "m")
ax1[1].plot(x_2, np.abs(v_cdm_2), label = rf"$k = {k_2.value:.2e}/\mathrm{{Mpc}}$", color = "b")
ax1[1].plot(x_1, np.abs(v_cdm_1), label = rf"$k = {k_1.value:.2e}/\mathrm{{Mpc}}$", color = "r")
ax1[1].plot(x_0, np.abs(v_cdm_0), label = rf"$k = {k_0.value:.2e}/\mathrm{{Mpc}}$", color = "g")

ax1[1].plot(x_5, np.abs(v_b_5), linestyle = "--", color = "orange")
ax1[1].plot(x_4, np.abs(v_b_4), linestyle = "--", color = "k")
ax1[1].plot(x_3, np.abs(v_b_3), linestyle = "--", color = "m")
ax1[1].plot(x_2, np.abs(v_b_2), linestyle = "--", color = "b")
ax1[1].plot(x_1, np.abs(v_b_1), linestyle = "--", color = "r")
ax1[1].plot(x_0, np.abs(v_b_0), linestyle = "--", color = "g")
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
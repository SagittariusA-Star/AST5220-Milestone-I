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

def Delta_CDM(delta, v, k, Hp):
    return delta - 3 * Hp / (const.c * k) * v


# Density parameters needed to plot background 
# color corresponding to domination era
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

eq_index  = np.argmin(np.abs(OmegaR - OmegaCDM + OmegaB))
k_horizon = 1 / eta_cosmo #(Hp_cosmo / const.c).to(1 / u.Mpc)
k_eq      = k_horizon[eq_index]
print("k_eq = : ", k_eq)

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
pertub_data_0 = np.loadtxt("perturbations_k5e-5.txt")
pertub_data_1 = np.loadtxt("perturbations_k2.85e-4.txt")
pertub_data_2 = np.loadtxt("perturbations_k1.623e-3.txt")
pertub_data_3 = np.loadtxt("perturbations_k9.244e-3.txt")
pertub_data_4 = np.loadtxt("perturbations_k5.266e-2.txt")
pertub_data_5 = np.loadtxt("perturbations_k3e-1.txt")

# Organizing data into individual arrays
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
Delta_CDM_0       = pertub_data_0[:, 11] 

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
Delta_CDM_1       = pertub_data_1[:, 11] 

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
Delta_CDM_2       = pertub_data_2[:, 11] 

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
Delta_CDM_3       = pertub_data_3[:, 11] 

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
Delta_CDM_4       = pertub_data_4[:, 11] 

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
Delta_CDM_5       = pertub_data_5[:, 11] 

# Plotting monopole and dipole moments of radiation perturbation as well as
# metric perturbations.
fig, ax = plt.subplots(2, 2 , figsize=[1.5 * 7.1014, 1.5 * 7.1014 / 1.618])

# Monopole moments
ax[0, 0].plot(x_5, Theta0_5, label=rf"$k = {k_5.value:.2e}/\mathrm{{Mpc}}$", color = "orange")
ax[0, 0].plot(x_4, Theta0_4, label=rf"$k = {k_4.value:.2e}/\mathrm{{Mpc}}$", color = "k")
ax[0, 0].plot(x_3, Theta0_3, label=rf"$k = {k_3.value:.2e}/\mathrm{{Mpc}}$", color = "m")
ax[0, 0].plot(x_2, Theta0_2, label=rf"$k = {k_2.value:.2e}/\mathrm{{Mpc}}$", color = "b")
ax[0, 0].plot(x_1, Theta0_1, label=rf"$k = {k_1.value:.2e}/\mathrm{{Mpc}}$", color = "r")
ax[0, 0].plot(x_0, Theta0_0, label=rf"$k = {k_0.value:.2e}/\mathrm{{Mpc}}$", color = "g")
#ax[0, 0].set_xlim(-12, 0)
ax[0, 0].set_ylabel(r"$\Theta_0$")
ax[0, 0].set_xlabel(r"$x = \log (a)$")
fig.legend(loc='center', bbox_to_anchor=(0.3, 0.3), fontsize = 14)

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

ax[0, 0].axvline(x_cosmo[np.argmin(np.abs(k_horizon - k_1))], color = "r", linestyle = ":")
ax[0, 0].axvline(x_cosmo[np.argmin(np.abs(k_horizon - k_2))], color = "b", linestyle = ":")
ax[0, 0].axvline(x_cosmo[np.argmin(np.abs(k_horizon - k_3))], color = "m", linestyle = ":")
ax[0, 0].axvline(x_cosmo[np.argmin(np.abs(k_horizon - k_4))], color = "k", linestyle = ":")
ax[0, 0].axvline(x_cosmo[np.argmin(np.abs(k_horizon - k_5))], color = "orange", linestyle = ":")
ax[0, 0].axvline(x_cosmo[np.argmin(np.abs(k_horizon - k_eq))], color = "r", linestyle = "-.")

# Dipole moments
ax[0, 1].plot(x_5, Theta1_5, label=rf"$k = {k_5.value:.2e}/\mathrm{{Mpc}}$", color = "orange")
ax[0, 1].plot(x_4, Theta1_4, label=rf"$k = {k_4.value:.2e}/\mathrm{{Mpc}}$", color = "k")
ax[0, 1].plot(x_3, Theta1_3, label=rf"$k = {k_3.value:.2e}/\mathrm{{Mpc}}$", color = "m")
ax[0, 1].plot(x_2, Theta1_2, label=rf"$k = {k_2.value:.2e}/\mathrm{{Mpc}}$", color = "b")
ax[0, 1].plot(x_1, Theta1_1, label=rf"$k = {k_1.value:.2e}/\mathrm{{Mpc}}$", color = "r")
ax[0, 1].plot(x_0, Theta1_0, label=rf"$k = {k_0.value:.2e}/\mathrm{{Mpc}}$", color = "g")
ax[0, 1].set_ylabel(r"$\Theta_1$ [c]")
ax[0, 1].set_xlabel(r"$x = \log (a)$")
#ax[0, 1].set_xlim(-12, 0)

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
ax[0, 1].axvline(x_cosmo[np.argmin(np.abs(k_horizon - k_1))], color = "r", linestyle = ":")
ax[0, 1].axvline(x_cosmo[np.argmin(np.abs(k_horizon - k_2))], color = "b", linestyle = ":")
ax[0, 1].axvline(x_cosmo[np.argmin(np.abs(k_horizon - k_3))], color = "m", linestyle = ":")
ax[0, 1].axvline(x_cosmo[np.argmin(np.abs(k_horizon - k_4))], color = "k", linestyle = ":")
ax[0, 1].axvline(x_cosmo[np.argmin(np.abs(k_horizon - k_5))], color = "orange", linestyle = ":")
ax[0, 1].axvline(x_cosmo[np.argmin(np.abs(k_horizon - k_eq))], color = "r", linestyle = "-.")

# Making space for legend of all plots
ax[1, 0].set_frame_on(False)
ax[1, 0].get_xaxis().set_visible(False)
ax[1, 0].get_yaxis().set_visible(False)

# Metric perturbations
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
#ax[1, 1].set_xlim(-12, 0)
ax[1, 1].set_ylabel(r"$\Phi$,  $\Psi$")
ax[1, 1].set_xlabel(r"$x = \log (a)$")

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
ax[1, 1].axvline(x_cosmo[np.argmin(np.abs(k_horizon - k_1))], color = "r", linestyle = ":")
ax[1, 1].axvline(x_cosmo[np.argmin(np.abs(k_horizon - k_2))], color = "b", linestyle = ":")
ax[1, 1].axvline(x_cosmo[np.argmin(np.abs(k_horizon - k_3))], color = "m", linestyle = ":")
ax[1, 1].axvline(x_cosmo[np.argmin(np.abs(k_horizon - k_4))], color = "k", linestyle = ":")
ax[1, 1].axvline(x_cosmo[np.argmin(np.abs(k_horizon - k_5))], color = "orange", linestyle = ":")
ax[1, 1].axvline(x_cosmo[np.argmin(np.abs(k_horizon - k_eq))], color = "r", linestyle = "-.")

fig.tight_layout()
fig.savefig("../doc/Figures/fig1.pdf", dpi=1000)

# Plotting density and velocity perturbations of matter components
fig1 = plt.figure(figsize=[1.5 * 7.1014, 1.5 * 7.1014 / 1.618])
ax10 = plt.subplot(221)
# CDM density perturbations
ax10.plot(x_5, np.abs(delta_cdm_5), label = fr"$k = {k_5.value:.2e}/\mathrm{{Mpc}}$", color = "orange")
ax10.plot(x_4, np.abs(delta_cdm_4), label = fr"$k = {k_4.value:.2e}/\mathrm{{Mpc}}$", color = "k")
ax10.plot(x_3, np.abs(delta_cdm_3), label = fr"$k = {k_3.value:.2e}/\mathrm{{Mpc}}$", color = "m")
ax10.plot(x_2, np.abs(delta_cdm_2), label = fr"$k = {k_2.value:.2e}/\mathrm{{Mpc}}$", color = "b")
ax10.plot(x_1, np.abs(delta_cdm_1), label = fr"$k = {k_1.value:.2e}/\mathrm{{Mpc}}$", color = "r")
ax10.plot(x_0, np.abs(delta_cdm_0), label = fr"$k = {k_0.value:.2e}/\mathrm{{Mpc}}$", color = "g")

# Baryon density perturbations
ax10.plot(x_5, np.abs(delta_b_5), linestyle = "--", color = "orange")
ax10.plot(x_4, np.abs(delta_b_4), linestyle = "--", color = "k")
ax10.plot(x_3, np.abs(delta_b_3), linestyle = "--", color = "m")
ax10.plot(x_2, np.abs(delta_b_2), linestyle = "--", color = "b")
ax10.plot(x_1, np.abs(delta_b_1), linestyle = "--", color = "r")
ax10.plot(x_0, np.abs(delta_b_0), linestyle = "--", color = "g")
#ax10.set_xlim(-12, 0)
#ax10.legend()
ax10.set_ylabel(r"$\delta_\mathrm{CDM}$, $|\delta_\mathrm{b}|$")
ax10.set_xlabel(r"$x = \log (a)$")
ax10.set_yscale("log")

ax10.axvspan(
    np.min(x_cosmo),
    x_cosmo[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    alpha=0.3,
    color="orange",
)
ax10.axvspan(
    x_cosmo[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    x_cosmo[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    alpha=0.3,
    color="cyan",
)
ax10.axvspan(
    x_cosmo[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    np.max(x_cosmo),
    alpha=0.3,
    color="purple",
)

ax10.axvspan(
    x_cosmo[recombo_index][0],
    x_cosmo[recombo_index][-1],
    alpha=0.4,
    color="red",
)

ax10.axvline(x_cosmo[np.argmin(np.abs(k_horizon - k_1))], color = "r", linestyle = ":")
ax10.axvline(x_cosmo[np.argmin(np.abs(k_horizon - k_2))], color = "b", linestyle = ":")
ax10.axvline(x_cosmo[np.argmin(np.abs(k_horizon - k_3))], color = "m", linestyle = ":")
ax10.axvline(x_cosmo[np.argmin(np.abs(k_horizon - k_4))], color = "k", linestyle = ":")
ax10.axvline(x_cosmo[np.argmin(np.abs(k_horizon - k_5))], color = "orange", linestyle = ":")
ax10.axvline(x_cosmo[np.argmin(np.abs(k_horizon - k_eq))], color = "r", linestyle = "-.")

# CDM velocity perturbations
ax11 = plt.subplot(222)
ax11.plot(x_5, np.abs(v_cdm_5), label = rf"$k = {k_5.value:.2e}/\mathrm{{Mpc}}$", color = "orange")
ax11.plot(x_4, np.abs(v_cdm_4), label = rf"$k = {k_4.value:.2e}/\mathrm{{Mpc}}$", color = "k")
ax11.plot(x_3, np.abs(v_cdm_3), label = rf"$k = {k_3.value:.2e}/\mathrm{{Mpc}}$", color = "m")
ax11.plot(x_2, np.abs(v_cdm_2), label = rf"$k = {k_2.value:.2e}/\mathrm{{Mpc}}$", color = "b")
ax11.plot(x_1, np.abs(v_cdm_1), label = rf"$k = {k_1.value:.2e}/\mathrm{{Mpc}}$", color = "r")
ax11.plot(x_0, np.abs(v_cdm_0), label = rf"$k = {k_0.value:.2e}/\mathrm{{Mpc}}$", color = "g")
# Baryon density perturbations
ax11.plot(x_5, np.abs(v_b_5), linestyle = "--", color = "orange")
ax11.plot(x_4, np.abs(v_b_4), linestyle = "--", color = "k")
ax11.plot(x_3, np.abs(v_b_3), linestyle = "--", color = "m")
ax11.plot(x_2, np.abs(v_b_2), linestyle = "--", color = "b")
ax11.plot(x_1, np.abs(v_b_1), linestyle = "--", color = "r")
ax11.plot(x_0, np.abs(v_b_0), linestyle = "--", color = "g")
#ax11.set_xlim(-12, 0)
#ax11.legend()
ax11.set_ylabel(r"$|v_\mathrm{CDM}|$, $|v_\mathrm{b}|$ [$c$]")
ax11.set_xlabel(r"$x = \log (a)$")
ax11.set_yscale("log")

ax11.axvspan(
    np.min(x_cosmo),
    x_cosmo[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    alpha=0.3,
    color="orange",
)
ax11.axvspan(
    x_cosmo[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    x_cosmo[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    alpha=0.3,
    color="cyan",
)
ax11.axvspan(
    x_cosmo[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    np.max(x_cosmo),
    alpha=0.3,
    color="purple",
)
ax11.axvspan(
    x_cosmo[recombo_index][0],
    x_cosmo[recombo_index][-1],
    alpha=0.4,
    color="red",
)
ax11.axvline(x_cosmo[np.argmin(np.abs(k_horizon - k_1))], color = "r", linestyle = ":")
ax11.axvline(x_cosmo[np.argmin(np.abs(k_horizon - k_2))], color = "b", linestyle = ":")
ax11.axvline(x_cosmo[np.argmin(np.abs(k_horizon - k_3))], color = "m", linestyle = ":")
ax11.axvline(x_cosmo[np.argmin(np.abs(k_horizon - k_4))], color = "k", linestyle = ":")
ax11.axvline(x_cosmo[np.argmin(np.abs(k_horizon - k_5))], color = "orange", linestyle = ":")
ax11.axvline(x_cosmo[np.argmin(np.abs(k_horizon - k_eq))], color = "r", linestyle = "-.")


# Co-moving dark matter density contrast
ax12 = plt.subplot(212)

ax12.plot(x_5, Delta_CDM_5, label = rf"$k = {k_5.value:.2e}/\mathrm{{Mpc}}$", color = "orange")
ax12.plot(x_4, Delta_CDM_4, label = rf"$k = {k_4.value:.2e}/\mathrm{{Mpc}}$", color = "k")
ax12.plot(x_3, Delta_CDM_3, label = rf"$k = {k_3.value:.2e}/\mathrm{{Mpc}}$", color = "m")
ax12.plot(x_2, Delta_CDM_2, label = rf"$k = {k_2.value:.2e}/\mathrm{{Mpc}}$", color = "b")
ax12.plot(x_1, Delta_CDM_1, label = rf"$k = {k_1.value:.2e}/\mathrm{{Mpc}}$", color = "r")
ax12.plot(x_0, Delta_CDM_0, label = rf"$k = {k_0.value:.2e}/\mathrm{{Mpc}}$", color = "g")
ax12.set_ylim(1e-6, 2e3)
ax12.legend()
ax12.set_ylabel(r"$\Delta_\mathrm{CDM}$")
ax12.set_xlabel(r"$x = \log (a)$")
ax12.set_yscale("log")

ax12.axvspan(
    np.min(x_cosmo),
    x_cosmo[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    alpha=0.3,
    color="orange",
)
ax12.axvspan(
    x_cosmo[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    x_cosmo[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    alpha=0.3,
    color="cyan",
)
ax12.axvspan(
    x_cosmo[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    np.max(x_cosmo),
    alpha=0.3,
    color="purple",
)
ax12.axvspan(
    x_cosmo[recombo_index][0],
    x_cosmo[recombo_index][-1],
    alpha=0.4,
    color="red",
)
ax12.axvline(x_cosmo[np.argmin(np.abs(k_horizon - k_1))], color = "r", linestyle = ":")
ax12.axvline(x_cosmo[np.argmin(np.abs(k_horizon - k_2))], color = "b", linestyle = ":")
ax12.axvline(x_cosmo[np.argmin(np.abs(k_horizon - k_3))], color = "m", linestyle = ":")
ax12.axvline(x_cosmo[np.argmin(np.abs(k_horizon - k_4))], color = "k", linestyle = ":")
ax12.axvline(x_cosmo[np.argmin(np.abs(k_horizon - k_5))], color = "orange", linestyle = ":")
ax12.axvline(x_cosmo[np.argmin(np.abs(k_horizon - k_eq))], color = "r", linestyle = "-.")

fig1.tight_layout()
fig1.savefig("../doc/Figures/fig2.pdf", dpi=1000)
plt.show()

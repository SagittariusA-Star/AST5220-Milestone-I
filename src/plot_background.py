import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from astropy import constants as const
from astropy import units as u

plt.style.use("bmh")


def redshift(x):
    """
    Redshift as a function of log-scale factor
    ----------------------
    Parameters:
        x: array or float
            Log-scale factor
    Returns:
        z: array of float
            Cosmological redshift
    ----------------------
    """
    return np.exp(-x) - 1


def log_scale_factor(z):
    """
    Log-Scale factor as a function of redshift
    ----------------------
    Parameters:
        z: array of float
            Cosmological redshift
    Returns:
        x: array or float
            Log-scale factor
    ----------------------
    """
    return -np.log(1 + z)


def eta_r(Hp):
    """
    Conformal time for radiation domminated era
    ----------------------
    Parameters:
        Hp: array of float
            Scaled Hubble parameter aH, 
            astropy units time
    Returns:
        eta : array or float
            Conformal time in astopy 
            units distance
    ----------------------
    """
    eta = const.c / Hp
    return eta


def eta_m(Hp, eta_matter, Hp_matter):
    """
    Conformal time for matter domminated era
    ----------------------
    Parameters:
        Hp: array of float
            Scaled Hubble parameter aH, 
            astropy units time
        eta_matter: float
            Conformal time in astropy units distance
            when matter density parameter Omega_m = 
            OmegaB + Omega_CDM is close to unity.
        Hp_matter: float
            Scaled Hubble parameter when matter domminates,
            i.e. Omega_m is close to unity.
    Returns:
        eta : array or float
            Conformal time in astopy 
            units distance
    ----------------------
    """
    eta = eta_matter + 2 * const.c * (1 / Hp - 1 / Hp_matter)
    return eta


def eta_L(Hp, eta_L, Hp_L):
    """
    Conformal time for dark energy domminated era
    ----------------------
    Parameters:
        Hp: array of float
            Scaled Hubble parameter aH, 
            astropy units time
        eta_L: float
            Conformal time in astropy units distance
            when matter density parameter OmegaLambda
            is close to unity.
        Hp_L: float
            Scaled Hubble parameter when dark energy domminates,
            i.e. Omega_Lambda is close to unity.
    Returns:
        eta : array or float
            Conformal time in astropy 
            units distance
    ----------------------
    """
    eta = eta_L + const.c * (1 / Hp_L - 1 / Hp)
    return eta


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
data = np.loadtxt("cosmology.txt")
x = data[:, 0]
eta_of_x = (data[:, 1] * u.m).to(u.Gpc)
Hp_of_x = (data[:, 2] / u.s).to(u.km / (u.s * u.Mpc))
OmegaB = data[:, 3]
OmegaCDM = data[:, 4]
OmegaLambda = data[:, 5]
OmegaR = data[:, 6]
Omega_sum = OmegaB + OmegaCDM + OmegaLambda + OmegaR
Omega_m = OmegaCDM + OmegaB

dHpdx_of_x = (data[:, 7] / u.s).to(u.km / (u.s * u.Mpc))
ddHpddx_of_x = (data[:, 8] / u.s).to(u.km / (u.s * u.Mpc))
print(np.max(ddHpddx_of_x / Hp_of_x))
eta_radiation = (eta_r(Hp_of_x)).to(u.Gpc)
eta_matter = (
    eta_m(
        Hp_of_x,
        eta_of_x[np.where(Omega_m == np.max(Omega_m))][-1],
        Hp_of_x[np.where(Omega_m == np.max(Omega_m))][-1],
    )
).to(u.Gpc)
eta_Lambda = (
    eta_L(
        Hp_of_x,
        eta_of_x[np.where(OmegaLambda == np.max(OmegaLambda))][-1],
        Hp_of_x[np.where(OmegaLambda == np.max(OmegaLambda))][-1],
    )
).to(u.Gpc)

# Plotting conformal time
fig, ax = plt.subplots(2, 2, figsize=[7.1014, 7.1014 / 1.618])
ax[0, 0].plot(x, eta_of_x, label=r"$\eta(x)$")
ax[0, 0].plot(x, eta_radiation, "r--", label=r"$\eta_r(x)$")
ax[0, 0].plot(x, eta_matter, "g:", label=r"$\eta_m(x)$")
ax[0, 0].plot(x, eta_Lambda, "m-.", label=r"$\eta_\Lambda(x)$")
ax[0, 0].legend()
ax[0, 0].set_xlabel(r"$x = \log (a)$")
ax[0, 0].set_ylabel(r"$\eta$ [Gpc]")
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


# Plotting Hubble parameter
ax[0, 1].semilogy(x, Hp_of_x, label=r"$aH(x)$")
ax[0, 1].legend()
ax[0, 1].set_xlabel(r"$x = \log (a)$")
ax[0, 1].set_ylabel(r"$aH$ [$\mathrm{km}\mathrm{s}^{-1} \mathrm{Mpc}^{-1}$]")
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

# Plotting Hubble parameter
ax[1, 0].semilogy(x, Hp_of_x / np.exp(x), label=r"$H(x)$")
ax[1, 0].scatter(0, 70, color="r")
ax[1, 0].legend()
ax[1, 0].set_xlabel(r"$x = \log (a)$")
ax[1, 0].set_ylabel(r"$H(x)$ [$\mathrm{km}\mathrm{s}^{-1} \mathrm{Mpc}^{-1}$]")
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
ax[1, 0].scatter(0, 70, color="r")


ax[1, 1].loglog(
    redshift(x[np.where(x < 0)]),
    Hp_of_x[np.where(x < 0)] / np.exp(x[np.where(x < 0)]),
    label=r"$H(a)$",
)
ax[1, 1].legend()
ax[1, 1].set_xlabel(r"$z$")

ax[1, 1].set_ylabel(r"$H(z)$ [$\mathrm{km}\mathrm{s}^{-1} \mathrm{Mpc}^{-1}$]")

# ax2 = ax[1, 1].twiny()

ax[1, 1].axvspan(
    np.max(redshift(x)),
    redshift(x)[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    alpha=0.3,
    color="orange",
)
ax[1, 1].axvspan(
    redshift(x)[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    redshift(x)[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    alpha=0.3,
    color="cyan",
)
ax[1, 1].axvspan(
    redshift(x)[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    np.min(redshift(x)),
    alpha=0.3,
    color="purple",
)


# ax2.set_xticks(np.linspace(1e-8, 1, 6))
# ax2.set_xlabel(r'$a$', color = "r")
# plt.tick_params(axis = "x", labelcolor = "r")
# ax2.set_xlim(1, 1e-8)
# ax2.grid(True, color = "r")
fig.tight_layout()
fig.savefig("../doc/Figures/Eta_&_H_of_x.pdf", dpi=1000)
plt.show()

fonts = {
    "font.family": "serif",
    "axes.labelsize": 15,
    "font.size": 13,
    "legend.fontsize": 15,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
}

plt.rcParams.update(fonts)


# Plotting density parameters
fig1, ax1 = plt.subplots(1, 1, figsize=[7.1014, 7.1014 / 1.618])

ax1.plot(x, OmegaB, label=r"$\Omega_b$")
ax1.plot(x, OmegaCDM, label=r"$\Omega_{CDM}$")
ax1.plot(x, OmegaLambda, label=r"$\Omega_\Lambda$")
ax1.plot(x, OmegaR, label=r"$\Omega_r$")
ax1.plot(x, Omega_sum, label=r"$\Sigma\Omega$")
ax1.legend()
ax1.set_xlabel(r"$x = \log (a)$")
ax1.set_ylabel(r"$\Omega$")
ax1.set_xlim(np.min(x), np.max(x))
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

plt.show()

fig1.tight_layout()
fig1.savefig("../doc/Figures/Omegas_of_x.pdf", dpi=1000)


fig2, ax2 = plt.subplots(1, 2, figsize=[7.1014, 7.1014 / 1.618])
ax2[0].plot(x, dHpdx_of_x / Hp_of_x, label=r"$\frac{\mathcal{H}'(x)}{\mathcal{H}}$")
ax2[0].set_xlabel(r"$x = \log (a)$")
ax2[0].set_ylabel(r"$\frac{\mathcal{H}'(x)}{\mathcal{H}}$")
ax2[0].axvspan(
    np.min(x),
    x[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    alpha=0.3,
    color="orange",
)
ax2[0].axvspan(
    x[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    x[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    alpha=0.3,
    color="cyan",
)
ax2[0].axvspan(
    x[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    np.max(x),
    alpha=0.3,
    color="purple",
)

ax2[0].axhline(0.5* (2 - 4), color = "r", label = r"$1 - \frac{b_r}{2}$")
ax2[0].axhline(0.5 * (2 - 3), linestyle = "--", color = "m", label = r"$1 - \frac{b_m}{2}$")
ax2[0].axhline(0.5 * 2, linestyle = "--", color = "g", label = r"$1 - \frac{b_\Lambda}{2}$")

ax2[1].plot(x, ddHpddx_of_x / Hp_of_x, label=r"$\frac{\mathcal{H}''(x)}{\mathcal{H}}$")
ax2[1].set_xlabel(r"$x = \log (a)$")
ax2[1].set_ylabel(r"$\frac{\mathcal{H}''(x)}{\mathcal{H}}$")
ax2[1].axvspan(
    np.min(x),
    x[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    alpha=0.3,
    color="orange",
)
ax2[1].axvspan(
    x[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0],
    x[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    alpha=0.3,
    color="cyan",
)
ax2[1].axvspan(
    x[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0],
    np.max(x),
    alpha=0.3,
    color="purple",
)


ax2[1].axhline(0.25* (2 - 4) ** 2, color = "r", label = r"$1 - \frac{b_r}{2}$")
ax2[1].axhline(0.25 * (2 - 3) ** 2, linestyle = "--", color = "m", label = r"$1 - \frac{b_m}{2}$")
ax2[1].axhline(0.25 * 2 ** 2, linestyle = "-.", color = "g", label = r"$1 - \frac{b_\Lambda}{2}$")
#ax2[1].set_yscale("log")

ax2[0].legend()
ax2[1].legend()
fig2.tight_layout()
fig2.savefig("../doc/Figures/derivatives.pdf", dpi=1000)

plt.show()

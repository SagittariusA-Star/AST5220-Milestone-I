import numpy as np 
import matplotlib.pyplot as plt
import matplotlib
from astropy import constants as const
from astropy import units as u

plt.style.use('bmh')




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
    return - np.log(1 + z) 


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
dHpdx_of_x = (data[:, 3] / u.s).to(u.km / (u.s * u.Mpc)) 
OmegaB      = data[:, 4]
OmegaCDM    = data[:, 5]
OmegaLambda = data[:, 6]
OmegaR      = data[:, 7]
Omega_sum = OmegaB + OmegaCDM + OmegaLambda + OmegaR

# Plotting conformal time
fig, ax = plt.subplots(2, 2, figsize=[7.1014, 7.1014 / 1.618])
ax[0, 0].plot(x, eta_of_x, label=r"$\eta(x)$")
ax[0, 0].legend()
ax[0, 0].set_xlabel(r"$x = \log (a)$")
ax[0, 0].set_ylabel(r"$\eta$ [Gpc]")

# Plotting Hubble parameter
ax[0, 1].semilogy(x, Hp_of_x, label=r"$aH(x)$")
ax[0, 1].legend()
ax[0, 1].set_xlabel(r"$x = \log (a)$")
ax[0, 1].set_ylabel(r"$aH$ [$\mathrm{km}\mathrm{s}^{-1} \mathrm{Mpc}^{-1}$]")

# Plotting Hubble parameter
ax[1, 0].semilogy(x, Hp_of_x / np.exp(x), label=r"$H(x)$")
ax[1, 0].legend()
ax[1, 0].set_xlabel(r"$x = \log (a)$")
ax[1, 0].set_ylabel(r"$H(x)$ [$\mathrm{km}\mathrm{s}^{-1} \mathrm{Mpc}^{-1}$]")

ax[1, 1].loglog(redshift(x[:-1]), Hp_of_x[:-1] / np.exp(x[:-1]), label=r"$H(a)$")
ax[1, 1].legend()
ax[1, 1].set_xlabel(r"$z$")

ax[1, 1].set_ylabel(r"$H(z)$ [$\mathrm{km}\mathrm{s}^{-1} \mathrm{Mpc}^{-1}$]")

ax2 = ax[1, 1].twiny()  
ax[1, 1].set_xlim(np.max(redshift(x[:-1])), np.min(redshift(x[:-1])))

ax2.set_xticks(np.linspace(1e-8, 1, 6))
ax2.set_xlabel(r'$a$', color = "r")
plt.tick_params(axis = "x", labelcolor = "r")  
ax2.set_xlim(1, 1e-8)
ax2.grid(True, color = "r")
fig.tight_layout()

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
ax1.plot(x, Omega_sum  , label=r"$\Sigma\Omega$")
ax1.legend()
ax1.set_xlabel(r"$x = \log (a)$")
ax1.set_ylabel(r"$\Omega$")
ax1.set_xlim(np.min(x), np.max(x))
ax1.axvspan(np.min(x), x[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0], alpha=0.3, color='orange')
ax1.axvspan(x[np.where(OmegaB + OmegaCDM >= OmegaLambda + OmegaR)][0], x[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0], alpha=0.3, color='cyan')
ax1.axvspan(x[np.where(OmegaLambda >= Omega_sum - OmegaLambda)][0], np.max(x), alpha=0.3, color='purple')


fig1.tight_layout()
#fig.savefig("../Figures/eta_of_x.pdf", dpi=1000)

plt.show()


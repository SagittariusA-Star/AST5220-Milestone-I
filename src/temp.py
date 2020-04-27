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
pertub_data_0 = np.loadtxt("perturbations_k0.1.txt")
x0 = pertub_data_0[:, 0]
S_j_ell = pertub_data_0[:, 12]

plt.plot(x0, S_j_ell)
plt.savefig("../doc/Figures/j_1200.pdf")
plt.show()

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


data = np.loadtxt("cells.txt")
ells = data[:, 0]
Cells = data[:, 1]
Theta_ell01 = data[:, 2]
Theta_ell001 = data[:, 3]
Theta_ell0001 = data[:, 4]


data1 = np.loadtxt("matter_power_spec.txt")
k = (data1[:, 0] / u.m).to(1 / u.Mpc)
P = (data1[:, 1] * u.m ** 3).to(u.Mpc ** 3)

fig = plt.figure(figsize=[1.5 * 7.1014, 1.5 * 7.1014 / 1.618])
ax = plt.subplot(111)
ax.plot(ells, Cells)
ax.set_yscale("log")
ax.set_xscale("log")

fig2 = plt.figure(figsize=[1.5 * 7.1014, 1.5 * 7.1014 / 1.618])
ax2 = plt.subplot(111)
ax2.plot(ells, Theta_ell01, label = r"$k = 0.1/\mathrm{Mpc}$")
ax2.plot(ells, Theta_ell001, label = r"$k = 0.01/\mathrm{Mpc}$")
ax2.plot(ells, Theta_ell0001, label = r"$k = 0.00/1\mathrm{Mpc}$")
ax2.set_yscale("symlog")
ax2.set_xscale("log")


fig1 = plt.figure(figsize=[1.5 * 7.1014, 1.5 * 7.1014 / 1.618])
ax1 = plt.subplot(111)
ax1.plot(k, P)
ax1.set_yscale("log")
ax1.set_xscale("log")
ax1.set_xlabel("k")
ax1.set_ylabel("P_M")

plt.show()

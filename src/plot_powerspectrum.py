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


data1 = np.loadtxt("matter_power_spec.txt")
k = (data1[:, 0] / u.m).to(1 / u.Mpc)
P = (data1[:, 1] * u.m ** 3).to(u.Mpc ** 3)
print(k)
fig = plt.figure(figsize=[1.5 * 7.1014, 1.5 * 7.1014 / 1.618])
ax = plt.subplot(111)
ax.plot(ells, Cells)
ax.set_yscale("symlog")


fig1 = plt.figure(figsize=[1.5 * 7.1014, 1.5 * 7.1014 / 1.618])
ax1 = plt.subplot(111)
ax1.plot(k, P)
ax1.set_yscale("log")
ax1.set_xscale("log")
ax1.set_xlabel("k")
ax1.set_ylabel("P_M")

plt.show()

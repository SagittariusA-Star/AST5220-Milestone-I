import numpy as np 
import matplotlib.pyplot as plt  


data = np.loadtxt("perturbations_k0.01.txt")
x = data[:, 0]
ST = data[:, -1]

plt.plot(x, ST, label = r"$S(x, k = 0.01 / Mpc)$")
plt.legend()
plt.xlabel("x")
plt.ylabel("S(x, k)")
plt.savefig("../doc/Figures/ST_k0.01.pdf", dpi=1000)

plt.show()
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import ellipk
from scipy.integrate import simpson


def dos_analytic(E, t=1.0):
    """Analytical DOS per unit cell (spin included) for t' = 0.
 
    E and t in the same units; returns rho in units of 1/[t]."""
    E = np.asarray(E, dtype=float)
    x = np.abs(E / t)
    rho = np.zeros_like(x)
 
    inner = x <= 1.0                       # -t <= E <= t
    outer = (x > 1.0) & (x <= 3.0)         # t < |E| <= 3t
 
    a = (1.0 + x) ** 2 - (x**2 - 1.0) ** 2 / 4.0
    b = 4.0 * x
 
    Z0 = np.where(inner, a, b)
    Z1 = np.where(inner, b, a)
 
    band = inner | outer
    with np.errstate(divide="ignore", invalid="ignore"):
        m = np.clip(Z1[band] / Z0[band], 0.0, 1.0)   # m = k^2
        rho[band] = (4.0 / np.pi**2) * (x[band] / abs(t)) / np.sqrt(Z0[band]) * ellipk(m)
    return rho


data = np.loadtxt('dos_nonorth.txt')
#for i in range(1, np.shape(data)[1]):
#    plt.plot(data[:,0], data[:,i], label='centered')
plt.plot(data[:,0], np.sum(data[:,1:], axis=1))
#plt.plot(data[:,0], dos_analytic(data[:,0], t=2.8), label='analytic')
data1 = np.loadtxt('dos_orth.txt')
plt.plot(data1[:,0], np.sum(data1[:,1:], axis=1))
plt.legend()
plt.show()
print(simpson(data[:,1], data[:,0]))
print(simpson(data1[:,1], data1[:,0]))


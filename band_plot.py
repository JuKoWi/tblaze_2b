import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('bands.txt')
plt.plot(data[:,3])
plt.plot(data[:,4])
plt.show()

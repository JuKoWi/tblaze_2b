import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('bands.txt')
for i in range(3, np.shape(data)[1]):
    plt.plot(data[:,i], '.', c='blue')
#for i in range(9,11):
#    plt.plot(data[:,i], '.', c='blue')
plt.show()

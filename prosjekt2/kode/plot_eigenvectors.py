import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import sys

data = np.loadtxt("../data/jacobi_eig2_10.dat")

n = np.arange(1,11,1)

eigenvalues = data[:,0]
eigenvectors = data[:,1:]

inds = eigenvalues.argsort()
sortedeigenvectors = eigenvectors[inds]

print(eigenvalues)
print(eigenvalues[inds])

for i in range(10):
    print(np.sum(eigenvectors[i]*eigenvectors[i]))
    plt.plot(n,eigenvectors[i])
plt.show()

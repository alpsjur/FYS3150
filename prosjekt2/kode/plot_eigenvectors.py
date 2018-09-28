import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import sys

def remove_file(filename):
    if os.path.exists(filename):
        os.system("rm " + filename)

def make_plot(ax, problem, n, omega_r,which_to_plot):

    os.system("./main.exe {} {} {}".format(n, problem, omega_r))
    #os.system("./main.exe {} {} 1 {}".format(n, problem, omega_r))

    filename = "../data/jacobi_eig{}_{}.dat".format(problem, n)
    data = np.loadtxt(filename)
    remove_file(filename)
    x = np.linspace(0,1,n)
    eigenvalues = data[:,0]
    eigenvectors = data[:,1:]
    inds = eigenvalues.argsort()
    sortedeigenvectors = eigenvectors[inds]

    for i in which_to_plot:
        ax.plot(x,eigenvectors[i]*n,label="n={}".format(n)) #ganger med n for å normalisere

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

for j in np.arange(10,101,10):
    make_plot(ax, 2, j, 5,[0])

ax.legend()
plt.show()

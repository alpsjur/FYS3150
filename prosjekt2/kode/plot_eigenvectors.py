import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import sys

def remove_file(filename):
    if os.path.exists(filename):
        os.system("rm " + filename)

def make_plot(ax, problem, n, omega_r, interact, which_to_plot):

    #os.system("./main.exe {} {} {}".format(n, problem, omega_r))
    os.system("./main.exe {} {} 1 {} {}".format(n, problem, omega_r, interact))

    filename = "../data/jacobi_eig{}_{}.dat".format(problem, n)
    data = np.loadtxt(filename)
    remove_file(filename)
    #x = np.linspace(1./n,1-1./n,n)
    x = np.linspace(0,1,n+2)
    eigenvalues = data[:,0]
    eigenvectors = data[:,1:]
    inds = eigenvalues.argsort()
    sortert = eigenvectors[inds]

    for i in which_to_plot:
        y = sortert[i]*(n+2)/np.sum(np.abs(sortert[i]))   #normaliserer egenvektoren
        y = np.insert(y, [0,n], [0,0])
        if problem in [1,2]:
            y = y*y
        ax.plot(x, y,label=r"n={}, $\omega_r = {}$, interact = {}".format(n, omega_r, interact))

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

'''
for j in np.arange(10,101,10):
    make_plot(ax, 0, j, 0.01 ,0 ,[2])
'''
omegalist = [0.01, 0.5, 1.0, 5.0]

make_plot(ax, 2, 200, 0.01, 0, [0])
make_plot(ax, 2, 200, 0.01, 1, [0])
ax.legend()

plt.show()

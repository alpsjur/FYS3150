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
<<<<<<< HEAD
    os.system("./main.exe {} {} 1 {} {}".format(n, problem, omega_r, interact))
=======
    os.system("./main.exe {} {} 1 {}".format(n, problem, omega_r))
>>>>>>> e9b9b03c304ec5cc94365d2259df0b948343be73

    filename = "../data/jacobi_eig{}_{}.dat".format(problem, n)
    data = np.loadtxt(filename)
    remove_file(filename)
<<<<<<< HEAD
    x = np.linspace(0,10,n)
=======
    x = np.linspace(1./n,1-1./n,n)
>>>>>>> e9b9b03c304ec5cc94365d2259df0b948343be73
    eigenvalues = data[:,0]
    eigenvectors = data[:,1:]
    inds = eigenvalues.argsort()
    sorted = eigenvectors[inds]

    if problem in [1,2]:
        normalised = sorted*sorted*n/10.
    else:
        normalised = sorted*sorted*n #FLAGG

    for i in which_to_plot:
<<<<<<< HEAD
        ax.plot(x,sortedeigenvectors[i]*sortedeigenvectors[i],label=r"n={}, $\omega_r = {}$, interact = {}".format(n, omega_r, interact))
=======
        ax.plot(x,normalised[i],label="n={}".format(n)) #ganger med n for å normalisere
>>>>>>> e9b9b03c304ec5cc94365d2259df0b948343be73

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

<<<<<<< HEAD
"""
for j in np.arange(10,201,10):
    make_plot(ax, 2, j, 0.01,[0])
"""
omegalist = [0.01, 0.5, 1.0, 5.0]
=======
for j in np.arange(10,101,10):
    make_plot(ax, 0, j, 1,[0])
>>>>>>> e9b9b03c304ec5cc94365d2259df0b948343be73

make_plot(ax, 2, 200, 0.01, 0, [0])
make_plot(ax, 2, 200, 0.01, 1, [0])
ax.legend()
plt.show()

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
    x = np.linspace(0,10,n+2)
    eigenvalues = data[:,0]
    eigenvectors = data[:,1:]
    inds = eigenvalues.argsort()

    sorted_eigenvectors = eigenvectors[inds]

    for i in which_to_plot:
        y = sorted_eigenvectors[i]
        y = np.insert(y, [0,n], [0,0])   #legger til grensene
        #normaliserer løsningen
        if problem in [1,2]:
            y *= y                       #kvadrerer for å få sansynlighetsfordelingen
            y = y*(n+2)/(np.sum(np.abs(y))*10)
        else:
            y = y*(n+2)/np.sum(np.abs(y))
        ax.plot(x, y,label=r"n={}, $\omega_r = {}$, interact = {}".format(n, omega_r, interact))

#plotting
sns.set()
sns.set_style("whitegrid")
sns.set_palette("husl")
plt.rc('text', usetex=True)
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

'''
for j in np.arange(10,101,10):
    make_plot(ax, 0, j, 0.01 ,0 ,[2])
'''
omegalist = [0.01, 0.5, 1.0, 5.0]
for omega in omegalist:
    make_plot(ax, 2, 200, omega, 0, [0])
    ax.set_xlabel(r'$\rho$',fontsize=14)
    ax.set_ylabel(r'$\Psi ^2$',fontsize=14)

ax.legend(fontsize=12)
plt.savefig("../figurer/quantom_dots.pdf")

plt.show()

import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import sys

def remove_file(filename):
    if os.path.exists(filename):
        os.system("rm " + filename)

def make_plot(ax, problem, n, omegalist, interact,which_to_plot):
    eigenvalues = np.zeros(len(omegalist))
    i = 0
    for omega_r in omegalist:
        os.system("./main.exe {} {} 1 {} {}".format(n, problem, omega_r, interact))

        filename = "../data/jacobi_eig{}_{}.dat".format(problem, n)
        data = np.loadtxt(filename)
        remove_file(filename)
        eigenvalues_temp = np.sort(data[:,0])
        eigenvalues[i] = eigenvalues_temp[which_to_plot]
        i += 1
    ax.plot(omegalist,eigenvalues,'o')
    ax.set_xlabel(r'$\omega_r$',fontsize=14)
    ax.set_ylabel(r'$\lambda_0$',fontsize=14)

omegalist = [0.01, 0.5, 1.0, 5.0]

#plotting
sns.set()
sns.set_style("whitegrid")
sns.set_palette("husl")
plt.rc('text', usetex=True)
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
make_plot(ax, 2, 200, omegalist, 0, 0)
plt.show()

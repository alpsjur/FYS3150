import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import sys

def make_file(problem, n, omega_r, interact=0):
    omegastr = omega_r[0] + omega_r[2] + omega_r[3]  # necessary for filename formatting
    filename = "../data/jacobi_eig_{}_{}_{}_{}.dat".format(problem, n, omegastr, interact)
    if os.path.exists(filename) == False:
        os.system("./main.exe {} {} 1 {} {}".format(n, problem, float(omega_r), interact))
    return filename

def plot_eigenvectors(ax, filename, n, problem, which_to_plot, color=None):
    rhomin = 0
    rhomax = 10

    data = np.loadtxt(filename)

    x = np.linspace(rhomin, rhomax, n+2)
    eigenvalues = data[:,0]
    eigenvectors = data[:,1:]
    inds = eigenvalues.argsort()

    sorted_eigenvalues = eigenvalues[inds]
    sorted_eigenvectors = eigenvectors[inds]

    for i in which_to_plot:
        y = sorted_eigenvectors[i]
        y = np.insert(y, [0,n], [0,0])   #legger til grensene
        #normaliserer løsningen
        if problem in [1,2]:
            y *= y                       #kvadrerer for å få sansynlighetsfordelingen
            y = y*(n+2)/(np.sum(np.abs(y))*(rhomax - rhomin))
        else:
            y = y*(n+2)/np.sum(np.abs(y))
        ax.plot(x, y, color=color)


def plot_eigenvalues(ax, filename, omega_r, which_to_plot, color=None):

    data = np.loadtxt(filename)

    eigenvalues_temp = np.sort(data[:,0])
    eigenvalue = eigenvalues_temp[which_to_plot]

    ax.scatter(float(omega_r), eigenvalue, color=color)

#plotting
sns.set()
sns.set_style("whitegrid")
sns.set_palette("husl")
plt.rc('text', usetex=True)

fig, ax = plt.subplots(2, 1)

n= 200
problem = 2
omegalist = ["0.01", "0.50", "1.00", "5.00"]
legends = []
for omega in omegalist:
    legends.append(r"$\omega_r = {}$".format(omega))
    filename_nointeract = make_file(problem, n, omega, interact=0)
    filename_interact = make_file(problem, n, omega, interact=1)

    plot_eigenvectors(ax[0], filename_nointeract, n, problem, [0], color='r')
    plot_eigenvalues(ax[1], filename_nointeract, omega, 0, color='r')

    plot_eigenvectors(ax[0], filename_interact, n, problem, [0], color='b')
    plot_eigenvalues(ax[1], filename_interact, omega, 0, color='b')

ax[0].set_xlabel(r'$\rho$', fontsize=14)
ax[0].set_ylabel(r'$\psi_0 ^2$', fontsize=14)
#ax[0].legend(legends, fontsize=14, frameon=False)
ax[0].legend(['no interaction', 'interaction'], fontsize=14, frameon=False)

ax[1].set_xlabel(r'$\omega_r$', fontsize=14)
ax[1].set_ylabel(r'$\lambda_0$', fontsize=14)

#plt.subplots_adjust(hspace=0.35)

#plt.savefig("../figurer/quantom_dots_nointer.pdf")

plt.show()

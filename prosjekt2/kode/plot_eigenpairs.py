import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import sys

def remove_file(filename):
    if os.path.exists(filename):
        os.system("rm " + filename)

def plot_eigenvectors(ax, problem, n, omega_r, interact, which_to_plot):
    rhomin = 0
    rhomax = 10

    #os.system("./main.exe {} {} {}".format(n, problem, omega_r))
    os.system("./main.exe {} {} 1 {} {}".format(n, problem, omega_r, interact))

    filename = "../data/jacobi_eig{}_{}.dat".format(problem, n)
    data = np.loadtxt(filename)
    remove_file(filename)
    #x = np.linspace(1./n,1-1./n,n)
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
        ax.plot(x, y,label=r"$\omega_r = {}$".format(omega_r))


def plot_eigenvalues(ax, problem, n, omegalist, interact, which_to_plot):
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

#plotting
sns.set()
sns.set_style("whitegrid")
sns.set_palette("husl")
plt.rc('text', usetex=True)

fig, ax = plt.subplots(2, 1)

omegalist = [0.01, 0.5, 1.0, 5.0]
for omega in omegalist:
    plot_eigenvectors(ax[0], 2, 200, omega, 0, [0])
plot_eigenvalues(ax[1], 2, 200, omegalist, 0, 0)

ax[0].set_xlabel(r'$\rho$', fontsize=14)
ax[0].set_ylabel(r'$\psi_0 ^2$', fontsize=14)
ax[0].legend(fontsize=14, frameon=False)

ax[1].set_xlabel(r'$\omega_r$', fontsize=14)
ax[1].set_ylabel(r'$\lambda_0$', fontsize=14)

plt.subplots_adjust(hspace=0.35)

plt.savefig("../figurer/quantom_dots.pdf")

plt.show()

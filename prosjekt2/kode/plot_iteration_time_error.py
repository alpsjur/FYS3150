import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import sys

def remove_file(filename):
    if os.path.exists(filename):
        os.system("rm " + filename)

#first command line argument sets the number of runns (N) for each value of n
N = int(sys.argv[1])
#second command line argument sets the maximum n
n_max = int(sys.argv[2])
#making a list of n values
n = np.arange(2,n_max+1)


#removing old data-logs
remove_file("../data/jacobi_log.dat")

#running the c++ programs
for i in n:
    for j in range(N):
        os.system("./main.exe {} {} {}".format(i, 0, 0))

#loading data from file
jacobi_data = np.loadtxt("../data/jacobi_log.dat")

jacobi_n = jacobi_data[:,0]
jacobi_iterations = jacobi_data[:,1]
jacobi_time = jacobi_data[:,2]
arma_time = jacobi_data[:,3]
max_error = jacobi_data[:,4]

#calculating average CPU-time and the standard deviation
jacobi_avgtime = np.zeros(n_max-1)
jacobi_std = np.zeros(n_max-1)

arma_avgtime = np.zeros(n_max-1)
arma_std = np.zeros(n_max-1)

for i in n:
    start = int((i-2)*N)
    stop = int((i-1)*N)

    tempj = jacobi_time[start:stop]
    tempa = arma_time[start:stop]

    jacobi_avgtime[int(i)-2] = np.mean(tempj)
    jacobi_std[int(i)-2] = np.std(tempj)

    arma_avgtime[int(i)-2] = np.mean(tempa)
    arma_std[int(i)-2] = np.std(tempa)

#plotting
sns.set()
sns.set_style("whitegrid")
sns.set_palette("husl")
plt.rc('text', usetex=True)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(np.log10(jacobi_n),np.log10(jacobi_iterations),
        zorder=10)
p = np.polyfit(np.log10(jacobi_n),np.log10(jacobi_iterations),1)
ax.plot(np.log10(n),np.polyval(p,np.log10(n)),
        linestyle='dashed',
        color='grey',
        zorder=5)
ax.text(1.5, 1.5, r'slope={:.2f}'.format(p[0]),
        fontsize=12,
        color='grey')
ax.set_xlabel(r'$\log_{10}$ n', fontsize=14)
ax.set_ylabel(r'$\log_{10}$ iterations', fontsize=14)
#plt.savefig("../figurer/n_vs_itterations.pdf")


fig2 = plt.figure()
ax2 = fig2.add_subplot(1, 1, 1)
ax2.scatter(np.log10(jacobi_n),np.log10(max_error))
ax2.set_xlabel(r'$\log_{10}$ n', fontsize=14)
ax2.set_ylabel(r'$\log_{10}$ $\epsilon_{max}$', fontsize=14)
#plt.savefig("../figurer/relative_error.pdf")


fig3 = plt.figure()
ax3 = fig3.add_subplot(1, 1, 1)
ax3.plot(np.log10(n), np.log10(jacobi_avgtime),label=r"Jacobis method")
ax3.plot(np.log10(n), np.log10(arma_avgtime),label=r"armadillo")
ax3.set_xlabel(r'$\log_{10}$ n', fontsize=14)
ax3.set_ylabel(r'$\log_{10}$ mean CPU-time [$\log_{10}$ms]', fontsize=14)
ax3.legend(fontsize=12)
#plt.savefig("../figurer/CPU_time.pdf")


plt.show()

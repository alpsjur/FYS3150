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

jacobi_data = np.loadtxt("../data/jacobi_log.dat")

n = jacobi_data[:,0]
jacobi_iterations = jacobi_data[:,1]
jacobi_time = jacobi_data[:,2]
arma_time = jacobi_data[:,3]
max_error = jacobi_data[:,4]

n_values = np.linspace(10, 100, 10)
for i in range(90, )
    arma_avgtime = np.mean(arma_time[-n_max:])
    arma_std = np.std(arma_time[-n_max:])

    jacobi_avgtime = np.mean(jacobi_time[-n_max:])
    jacobi_std = np.std(jacobi_time[-n_max:])


print("Mean CPU times and standard deviation")
print("Jacobi time = {}, Jacobi std = {}".format(jacobi_avgtime, jacobi_std))
print("Armadillo time = {}, Armadillo std = {}".format(arma_avgtime, arma_std))
print("Time ratio = {}".format(jacobi_avgtime/arma_avgtime))


sns.set()
sns.set_style("whitegrid")
sns.set_palette("husl")
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(np.log10(n),np.log10(jacobi_iterations))
ax.set_xlabel(r'$\log_{10}$ n', fontsize=14)
ax.set_ylabel(r'$\log_{10}$ iterations', fontsize=14)
#plt.savefig("../figurer/n_vs_itterations.pdf")


fig2 = plt.figure()
ax2 = fig2.add_subplot(1, 1, 1)
ax2.plot(np.log10(n),np.log10(max_error))
ax2.set_xlabel(r'$\log_{10}$ n', fontsize=14)
ax2.set_ylabel(r'$\log_{10}$ $\epsilon$', fontsize=14)
#plt.savefig("../figurer/relative_error.pdf")


plt.show()

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
        os.system("./jacobis_method.exe {}".format(i))

jacobi_data = np.loadtxt("../data/jacobi_log.dat")

jacobi_n = jacobi_data[:,0]
jacobi_itterations = jacobi_data[:,1]
jacobi_time = jacobi_data[:,2]
arma_time = jacobi_data[:,3]

sns.set()
sns.set_style("whitegrid")
sns.set_palette("husl")
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(np.log10(jacobi_n),np.log10(jacobi_itterations))
ax.set_xlabel(r'$\log_{10}$ n', fontsize=14)
ax.set_ylabel(r'$\log_{10}$ itterations', fontsize=14)
plt.savefig("../figurer/n_vs_itterations.pdf")
plt.show()

time_ratio = jacobi_time/arma_time

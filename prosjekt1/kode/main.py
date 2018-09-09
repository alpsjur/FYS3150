'''
main python script for running all the c++ programs, plotting the error and
printing the time ratio
'''
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
#second command line argument sets the maximum power of n
n_max = int(sys.argv[2])
#making a list of n power values
n = np.arange(1,n_max+1)

#ask the user if the LU decomposition should the run.
ans = input("Run LU decomposition? [y/n] ")
if ans == 'y':
    LU = True
else:
    LU = False

#removing old data-logs
remove_file("../data/general_matrix_time_log.dat")
remove_file("../data/special_matrix_time_log.dat")
remove_file("../data/LU_time_log.dat")
remove_file("../data/max_error_log.dat")

#running the c++ programs
for i in n:
    for j in range(N):
        os.system("./general_matrix_solver.exe {}".format(i))
        os.system("./special_matrix_solver.exe {}".format(i))
        if LU:
            os.system("./LU_decomp.exe {}".format(i))

#load the generated data from the datalogs
general_matrix_time = np.loadtxt("../data/general_matrix_time_log.dat")
special_matrix_time = np.loadtxt("../data/special_matrix_time_log.dat")
if LU:
    LU_time = np.loadtxt("../data/LU_time_log.dat")
max_error = np.loadtxt("../data/max_error_log.dat")

#plotting maximum error plot
sns.set()
sns.set_style("whitegrid")
sns.set_palette("husl")
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(max_error[:,0],max_error[:,1])
ax.set_xlabel(r'$\log_{10}$ n', fontsize=14)
ax.set_ylabel(r'$\epsilon$', fontsize=14)
plt.savefig("../figurer/error.pdf")
plt.show()

#calculation the ratio between the CPU time of the special and general algorithm,
#and the LU decomposition and special algorithm
time_diff = np.zeros((n_max,3))
for i in n:
    temp1 = general_matrix_time[(i-1)*N:i*N,1]/special_matrix_time[(i-1)*N:i*N,1]
    if LU:
        temp2 = LU_time[(i-1)*N:i*N,1]/special_matrix_time[(i-1)*N:i*N,1]
    else:
        temp2 = np.nan
    time_diff[i-1] = np.array([i,np.mean(temp1), np.mean(temp2)])

#printing the time ratio 
print(time_diff)

import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import sys

def remove_file(filename):
    if os.path.exists(filename):
        os.system("rm " + filename)

N = int(sys[1])
n_max = int(sys[2])
n = np.arange(1,n_max+1)

ans = input("Run LU decomposition? [y/n]")
if ans == 'y':
    LU = True
else:
    LU = False

remove_file("../data/general_matrix_time_log.dat")
remove_file("../data/special_matrix_time_log.dat")
remove_file("../data/LU_time_log.dat")
remove_file("../data/max_error_log.dat")

for i in n:
    for j in range(N):
        os.system("./general_matrix_solver.exe {}".format(i))
        os.system("./special_matrix_solver.exe {}".format(i))
        if LU:
            os.system("./LU_decomp.exe {}".format(i))


general_matrix_time = np.loadtxt("../data/general_matrix_time_log.dat")
special_matrix_time = np.loadtxt("../data/special_matrix_time_log.dat")
if LU:
    LU_time = np.loadtxt("../data/LU_time_log.dat")
max_error = np.loadtxt("../data/max_error_log.dat")

#plt.plot(max_error[:,0],max_error[:,1])
#plt.show()

#Regner ut forholdet mellom CPU-tid for general og special
time_diff = np.zeros((n_max,3))
for i in n:
    temp1 = general_matrix_time[(i-1)*N:i*N,1]/special_matrix_time[(i-1)*N:i*N,1]
    if LU:
        temp2 = general_matrix_time[(i-1)*N:i*N,1]/LU_time[(i-1)*N:i*N,1]
    else:
        temp2 = np.nan
    time_diff[i-1] = np.array([i,np.mean(temp1), np.mean(temp2)])

print(time_diff)

import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

'''
os.system("rm time.log")

N = 100
n = 4

for i in range(N):
    os.system("./general_matrix_solver.exe {} >> time.log".format(n))
'''

N = 10
n_max = 7
n = np.arange(1,n_max+1)


os.system("rm ../data/general_matrix_time_log.dat")
os.system("rm ../data/special_matrix_time_log.dat")
os.system("rm ../data/max_error_log.dat")

for i in n:
    for j in range(N):
        os.system("./general_matrix_solver.exe {}".format(i))
        os.system("./special_matrix_solver.exe {}".format(i))


general_matrix_time = np.loadtxt("../data/general_matrix_time_log.dat")
special_matrix_time = np.loadtxt("../data/special_matrix_time_log.dat")
max_error = np.loadtxt("../data/max_error_log.dat")

plt.plot(max_error[:,0],max_error[:,1])
plt.show()

#Regner ut relativ forskjell i CPU-tid
time_diff = np.zeros((n_max,2))
for i in n:
    temp = (general_matrix_time[(i-1)*N:i*N,1] - \
     special_matrix_time[(i-1)*N:i*N,1])/special_matrix_time[(i-1)*N:i*N,1]
    time_diff[i-1] = np.array([i,np.mean(temp)])

print(time_diff)

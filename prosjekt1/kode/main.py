import os

os.system("rm time.log")

N = 100
n = 4

for i in range(N):
    os.system("./general_matrix_solver.exe {} >> time.log".format(n))

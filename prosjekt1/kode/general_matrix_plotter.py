import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

def u(x):
    f = 1 - (1 - np.exp(-10))*x - np.exp(-10*x)
    return f


datalist = [10, 100, 1000]


sns.set()
sns.set_style("whitegrid")
sns.set_palette("husl")
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

legends = []
j = 0
for i in datalist:
    data = np.loadtxt("../data/general_matrix{}.dat".format(i))
    ax.plot(data[:, 0], data[:,1])
    legends.append('n = {}'.format(i))
    j += 1

ax.plot(data[:,0], u(data[:,0]))

legends.append('analytic')

ax.set_ylabel('u(x)', fontsize=14)
ax.set_xlabel("x", fontsize=14)
ax.legend(legends, loc="best", fontsize=12)


plt.savefig("../figurer/general_matrix_compare.pdf")

plt.show()

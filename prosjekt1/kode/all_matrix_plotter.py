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

fig, ax = plt.subplots(3, 1, sharex='col', sharey='row')

dat = ["general_matrix","special_matrix","LU_decomp"]
name = ["general algorithm", "special algorithm", "LU decomposition"]
for j in range(3):
    legends = []
    for i in datalist:
        data = np.loadtxt("../data/" + dat[j] + "{}.dat".format(i))
        ax[j].plot(data[:,0], data[:,1],label='n = {}'.format(i))
    ax[j].plot(data[:,0], u(data[:,0]),label='analytic')
    ax[j].text(0.3, 0.15, name[j], horizontalalignment='center', \
    verticalalignment='center', color="grey")


ax[0].legend(loc='upper center', bbox_to_anchor=(0.5, 1.3),
          ncol=4)



fig.text(0.5, 0.04, 'x', fontsize=14, ha='center')
fig.text(0.03, 0.5, 'u(x)', fontsize=14, va='center', rotation='vertical')


plt.savefig("../figurer/all_matrix_compare.pdf")

plt.show()

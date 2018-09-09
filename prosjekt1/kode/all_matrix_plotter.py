'''
Script for plotting the solutions for n=10,100,1000 coputed with the general,
specialized and LU decomposition algorithm.
'''
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

#analytic solution
def u(x):
    f = 1 - (1 - np.exp(-10))*x - np.exp(-10*x)
    return f

#values for n
datalist = [10, 100, 1000]

#setting plot style using seaborn
sns.set()
sns.set_style("whitegrid")
sns.set_palette("husl")

#figure with three subplots
fig, ax = plt.subplots(3, 1, sharex='col', sharey='row')

#list of file names and algorithm names
filename = ["general_matrix","special_matrix","LU_decomp"]
algorithm = ["general algorithm", "special algorithm", "LU decomposition"]

#plotting
for j in range(3):
    legends = []
    for i in datalist:
        data = np.loadtxt("../data/" + filename[j] + "{}.dat".format(i))
        ax[j].plot(data[:,0], data[:,1],label='n = {}'.format(i))
    #ax[j].plot(data[:,0], u(data[:,0]),label='analytic')
    ax[j].text(0.3, 0.15, algorithm[j], horizontalalignment='center', \
    verticalalignment='center', color="grey")

#adding legends above the first plot
ax[0].legend(loc='upper center', bbox_to_anchor=(0.5, 1.3),
          ncol=4,fontsize=12)

#adding one common x and y label
fig.text(0.5, 0.04, 'x', fontsize=14, ha='center')
fig.text(0.03, 0.5, 'u(x)', fontsize=14, va='center', rotation='vertical')


plt.savefig("../figurer/all_matrix_compare.pdf")

plt.show()

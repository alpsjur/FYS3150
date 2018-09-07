import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset

def u(x):
    f = 1 - (1 - np.exp(-10))*x - np.exp(-10*x)
    return f


datalist = [10, 100, 1000]


sns.set()
sns.set_style("whitegrid")
sns.set_palette("husl")
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
axins = zoomed_inset_axes(ax, 200, loc=8)

legends = []
j = 0
for i in datalist:
    data = np.loadtxt("../data/general_matrix{}.dat".format(i))
    ax.plot(data[:, 0], data[:,1])
    axins.plot(data[:, 0], data[:,1])
    legends.append('n = {}'.format(i))
    j += 1

ax.plot(data[:,0], u(data[:,0]))
axins.plot(data[:,0], u(data[:,0]))

legends.append('analytic')

ax.set_ylabel('v(x)', fontsize=14)
ax.set_xlabel("x", fontsize=14)
ax.legend(legends, loc=1, fontsize=12)

x1, x2, y1, y2 = 0.228, 0.229, 0.669, 0.67 # specify the limits
axins.set_xlim(x1, x2) # apply the x-limits
axins.set_ylim(y1, y2) # apply the y-limits
#plt.yticks(visible=False)
plt.xticks(visible=False)
mark_inset(ax, axins, loc1=1, loc2=2, fc="none", ec="0.5")
plt.savefig("../figurer/general_matrix_compare.pdf")

plt.show()

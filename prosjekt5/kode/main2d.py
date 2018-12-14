import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import cm
import seaborn as sns


sns.set()
sns.set_style("whitegrid")
sns.set_palette("husl")

datadir = "../data/"
figdir = "../figurer/"


data_centered = np.fromfile(datadir + "psi_bounded_sine_centered_2d.bin")
data_centered = np.reshape(data_centered, ((1500,41,41)))
#data_centered = np.fromfile(datadir + "psi_bounded_sine_centered_2d.bin")

fig, axes = plt.subplots(ncols=2, nrows=2, sharex=True, sharey=True)
fig.add_subplot(111, frameon=False)

plt.tick_params(labelcolor="none", top=False, bottom=False, left=False, right=False)
plt.grid(False)
plt.xlabel("Zonal extent", fontsize=14)
plt.ylabel("Meridional extent", fontsize=14)

x = np.linspace(0, 1, 41)
y = np.linspace(0, 1, 41)
x,y = np.meshgrid(x,y)

times = [50, 100, 150]
labels = ["t = 0"]
#ax[0,0].contourf(x, y, data_forward[0])
#ax[0,1].contourf(x, y, data_forward[499])
#ax[1,0].contourf(x, y, data_forward[999])
#ax[1,1].contourf(x, y, data_forward[1499])
#ax[1].contourf(x, y, data_centered[0])


for ax, tindex in zip(axes.flat, [0, 499, 999, 1499]):
    im = ax.contourf(x, y, data_centered[tindex], cmap=cm.Greens_r)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(im, cax=cbar_ax)
#plt.savefig(figdir + "bounded_sine_centered_2d.pdf")

plt.show()

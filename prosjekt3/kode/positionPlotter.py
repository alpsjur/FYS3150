from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import glob
import seaborn as sns


def plotPlanet(ax, filename, d3=False):
    data = np.loadtxt(filename)
    name = filename.split("/")[-1][:-4]
    x = data[:,0]; y = data[:,1]; z = data[:,2];
    if d3:
        ax.plot3D(x,y,z,label=name)
    else:
        ax.plot(x,y,label=name)

def plotSystem(ax, path, d3=False):
    files = glob.glob(path+"/*.dat")
    for filename in files:
        plotPlanet(ax, filename, d3)


directory = "../figurer/"
if not os.path.exists(directory):
    os.makedirs(directory)

sns.set()
sns.set_style("white")
sns.set_palette("husl")
#plt.rc('text', usetex=True)

fig = plt.figure()
ax = plt.axes(projection='3d')
#ax = fig.add_subplot(1,1,1)

plotSystem(ax, "../data/solarsystem")

ax.legend()
ax.grid(False)
plt.axis('equal')
#plt.axis('off')
plt.savefig(directory + "fig.pdf")
plt.show()

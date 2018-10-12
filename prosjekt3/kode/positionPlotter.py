from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np
import sys
import glob
import seaborn as sns

sns.set()
sns.set_style("white")
sns.set_palette("husl")
#plt.rc('text', usetex=True)

def plotPlanet(ax, filename):
    data = np.loadtxt(filename)
    name = filename.split("/")[-1][:-4]
    x = data[:,0]; y = data[:,1]; z = data[:,2];
    #ax.plot3D(x,y,z,label=name)
    ax.plot(x,y,label=name)

def plotSystem(ax, path):
    files = glob.glob(path+"/*.dat")
    for filename in files:
        plotPlanet(ax, filename)


fig = plt.figure()
#ax = plt.axes(projection='3d')
ax = fig.add_subplot(1,1,1)

<<<<<<< HEAD
plotPlanet(ax, "../data/solarsystem/Sun.dat")
plotPlanet(ax, "../data/solarsystem/Mercury.dat")
plotPlanet(ax, "../data/solarsystem/Venus.dat")
plotPlanet(ax, "../data/solarsystem/Earth.dat")
plotPlanet(ax, "../data/solarsystem/Mars.dat")
plotPlanet(ax, "../data/solarsystem/Jupiter.dat")
plotPlanet(ax, "../data/solarsystem/Saturn.dat")
plotPlanet(ax, "../data/solarsystem/Uranus.dat")
plotPlanet(ax, "../data/solarsystem/Neptune.dat")
=======

plotSystem(ax, "../data/sun_earth_jupiter/jupiter_mass_1000")
>>>>>>> 07f3b98035ff0d97acfee72655ea3463891b68f7

ax.legend()
ax.grid(False)
plt.axis('equal')
#plt.axis('off')

plt.show()

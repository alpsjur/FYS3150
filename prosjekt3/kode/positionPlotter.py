from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

sns.set()
sns.set_style("white")
sns.set_palette("husl")
plt.rc('text', usetex=True)

def plotPlanet(ax, filename):
    data = np.loadtxt(filename)
    name = filename[8:-4]
    x = data[:,0]; y = data[:,1]; z = data[:,2];
    #ax.plot3D(x,y,z,label=name)
    ax.plot(x,y,label=name)

fig = plt.figure()
#ax = plt.axes(projection='3d')
ax = fig.add_subplot(1,1,1)

plotPlanet(ax, "../data/Earth.dat")
plotPlanet(ax, "../data/Jupiter.dat")
plotPlanet(ax, "../data/Sun.dat")
plotPlanet(ax, "../data/Mars.dat")
plotPlanet(ax, "../data/Mercury.dat")
plotPlanet(ax, "../data/Neptune.dat")
plotPlanet(ax, "../data/Saturn.dat")
plotPlanet(ax, "../data/Uranus.dat")
plotPlanet(ax, "../data/Venus.dat")


ax.legend()
ax.grid(False)
plt.axis('equal')
#plt.axis('off')

plt.show()

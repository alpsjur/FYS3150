from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

sns.set()
sns.set_style("white")
sns.set_palette("husl")
plt.rc('text', usetex=True)

def plotPlanet(ax, filename):
    data = np.loadtext(filename)
    name = filename[:-4]
    x = data[0]; y = data[1]; z = data[2];
    ax.plot3D(x,y,z,label=name)

fig = plt.figure()
ax = plt.axes(projection='3d')


ax.legend()
ax.grid(False)
#plt.axis('off')

plt.show()

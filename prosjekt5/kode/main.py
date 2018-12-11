import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import seaborn as sns


sns.set()
#sns.set_style("whitegrid")
sns.set_palette("husl")

datadir = "../data/"
figdir = "../figurer/"
data = np.loadtxt(datadir + "psi_bounded.dat")
x = np.linspace(0, 1, len(data[0]))
t = np.linspace(0, 150, len(data[:, 0]))

ymax = np.max(data)
ymin = np.min(data)
addedspace = (ymax-ymin)/10

fig, ax = plt.subplots()
#ax.set_xlim(-0.1,1.1)
#ax.set_ylim(ymin-addedspace,ymax+addedspace)

#line, = ax.plot(x, data[0])

def init():
    return line,

def animate(i):
    line.set_ydata(data[i])
    return line,

"""
anim = animation.FuncAnimation(fig, animate,
                                frames=len(data),
                                interval=2,
                                blit=True
                                )
"""
c = ax.contourf(x, t, data)
ax.set_xlabel("Spatial extent")
ax.set_ylabel("Time")
fig.colorbar(c, label="Amplitude")
plt.savefig(figdir + "hovmuller_boundedgaussian.pdf")
plt.show()

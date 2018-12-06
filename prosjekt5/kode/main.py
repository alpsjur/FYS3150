import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import seaborn as sns


sns.set()
sns.set_style("whitegrid")
sns.set_palette("husl")


datadir = "../data/"
data = np.loadtxt(datadir + "psi.dat")
x = np.linspace(0, 1, len(data[0]))

fig, ax = plt.subplots()

line, = ax.plot(x, data[0])

def init():
    ax.set_xlim(-0.1, 1.1)
    #ax.set_ylim(-0.15, 0.15)
    return line,

def animate(i):
    line.set_ydata(data[i+1])
    return line,


anim = animation.FuncAnimation(fig, animate,
                                init_func=init,
                                frames=200,
                                interval=2,
                                blit=True
                                )
plt.show()

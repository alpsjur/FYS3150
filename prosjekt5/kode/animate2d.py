import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import seaborn as sns


sns.set()
sns.set_style("whitegrid")
sns.set_palette("husl")

datadir = "../data/"
file = datadir + "psi_periodic_gaussian_centered_2d.bin"
data = np.fromfile(file)
data = np.reshape(data,((1500,41,41)))
'''

x = np.linspace(0, 1, 41)
y = np.linspace(0, 1, 41)
x,y = np.meshgrid(x,y)

fig = plt.figure()
ax = plt.axes(xlim=(0, 1), ylim=(0, 1))
plt.xlabel(r'x')
plt.ylabel(r'y')

# animation function
def animate(i):
    z = data[i,:,:]
    #ax.clear()
    cont = ax.contourf(x, y, z, 25)
    return cont

anim = animation.FuncAnimation(fig, animate, frames=1500, interval = 2, blit=True)

plt.show()
'''
data = data[:,20,:]
x = np.linspace(0, 1, len(data[0]))

ymax = np.max(data)
ymin = np.min(data)
addedspace = (ymax-ymin)/10

fig, ax = plt.subplots()
ax.set_xlim(-0.1,1.1)
ax.set_ylim(ymin-addedspace,ymax+addedspace)

line, = ax.plot(x, data[0])

def init():
    return line,

def animate(i):
    line.set_ydata(data[i])
    return line,


anim = animation.FuncAnimation(fig, animate,
                                frames=len(data),
                                interval= 2,
                                blit=True
                                )
plt.show()

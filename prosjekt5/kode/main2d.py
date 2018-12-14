import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import seaborn as sns


sns.set()
sns.set_style("whitegrid")
sns.set_palette("husl")

datadir = "../data/"
figdir = "../figurer/"


data_forward = np.fromfile(datadir + "psi_periodic_gaussian_centered_2d.bin")
data_forward = np.reshape(data_forward,((1500,41,41)))
#data_centered = np.fromfile(datadir + "psi_bounded_sine_centered_2d.bin")

fig, ax = plt.subplots(2, 2, sharex=True)
fig.add_subplot(111, frameon=False)

plt.tick_params(labelcolor="none", top=False, bottom=False, left=False, right=False)
plt.grid(False)
plt.xlabel("Spatial extent", fontsize=14)
plt.ylabel("Amplitude", fontsize=14)

x = np.linspace(0, 1, 41)
y = np.linspace(0, 1, 41)
x,y = np.meshgrid(x,y)

times = [50, 100, 150]
labels = ["t = 0"]
ax[0,0].contourf(x, y, data_forward[0])
ax[0,1].contourf(x, y, data_forward[499])
ax[1,0].contourf(x, y, data_forward[999])
ax[1,1].contourf(x, y, data_forward[1499])
#ax[1].contourf(x, y, data_centered[0])

plt.show()

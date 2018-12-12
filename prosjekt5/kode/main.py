import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import seaborn as sns


sns.set()
sns.set_style("whitegrid")
sns.set_palette("husl")

datadir = "../data/"
figdir = "../figurer/"


data_forward = np.loadtxt(datadir + "psi_periodic_sine_forward.dat")
data_centered = np.loadtxt(datadir + "psi_periodic_sine_centered.dat")

x = np.linspace(0, 1, len(data_forward[0]))
t = np.linspace(0, 150, len(data_forward[:, 0]))
dt = t[-1]/len(data_forward[:, 0])

fig, ax = plt.subplots(2, 1, sharex=True)
fig.add_subplot(111, frameon=False)

plt.tick_params(labelcolor="none", top=False, bottom=False, left=False, right=False)
plt.grid(False)
plt.xlabel("Spatial extent", fontsize=14)
plt.ylabel("Amplitude", fontsize=14)

times = [50, 100, 150]
labels = ["t = 0"]
ax[0].plot(x, data_forward[0])
ax[1].plot(x, data_centered[0])
for time in times:
    labels.append("t = {}".format(time))
    timeIndex = int(time/dt) - 1
    ax[0].plot(x, data_forward[timeIndex])
    ax[1].plot(x, data_centered[timeIndex])
fig.legend(labels, frameon=False, ncol=4, bbox_to_anchor=(1.0, 0.96), fontsize=14)
#plt.savefig(figdir + "compare_forward_centered.pdf")

# generating hov muller diagram
fig, ax = plt.subplots()
c = ax.contourf(x, t, data_centered)
ax.set_xlabel("Spatial extent", fontsize=14)
ax.set_ylabel("Time", fontsize=14)
fig.colorbar(c, label="Amplitude")
#plt.savefig(figdir + "hovmuller_periodicgaussian.pdf")
plt.show()

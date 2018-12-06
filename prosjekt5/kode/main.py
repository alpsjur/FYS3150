import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


sns.set()
sns.set_style("whitegrid")
sns.set_palette("husl")


datadir = "../data/"
data = np.loadtxt(datadir + "psi.dat")
print(data)


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

for i in [33, 66, 99]:
    ax.plot(data[i])
plt.show()

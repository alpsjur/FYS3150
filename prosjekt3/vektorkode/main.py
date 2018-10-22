import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def remove_file(filedir):
    if os.path.exists(filedir):
        os.system("rm " + filedir)

def create_directory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def emptyax(fig):
    ax = fig.add_subplot(1, 1, 1)
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
    return ax


def run_maincpp(scenario, endtime, dt, *args, n=1):
    if scenario == 2:
        velocityScale = args[0]
        beta = args[1]
        for i in range(n):
            os.system("./main.exe {} {} {} {} {}".format(scenario, endtime, dt, velocityScale, beta))
    else:
        for i in range(n):
            os.system("./main.exe {} {} {}".format(scenario, endtime, dt))

def plotPlanet(ax, filename, d3=False):
    data = np.loadtxt(filename)
    name = filename.split("/")[-1][:-4]
    x = data[:, 0]; y = data[:, 1]; z = data[:, 2]
    if d3:
        ax.plot3D(x, y, z, label=name)
    else:
        ax.plot(x, y, label=name)

def plotSystem(ax, path, d3=False):
    files = glob.glob(path+"/*.dat")
    for filename in files:
        plotPlanet(ax, filename, d3)


if __name__ == "__main__":
    # sammenligne Euler og Verlet
    path = "../data/euler_vs_verlet"
    fig = plt.figure()
    ax = fig.add_subplot(3, 1, 1)

    run_maincpp(1, 1, 0.001)
    plotSystem(ax, path)
    run_maincpp(1, 5, 0.001)
    ax = fig.add_subplot(3, 1, 2)
    plotSystem(ax, path)
    run_maincpp(1, 10, 0.001)
    ax = fig.add_subplot(3, 1, 3)
    plotSystem(ax, path)

    ax = emptyax(fig)
    ax.set_xlabel("x [AU]")
    ax.set_ylabel("y [AU]")
    plt.show()

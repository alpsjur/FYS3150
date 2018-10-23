import os
import sys
import glob
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
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

def plotPlanet(ax, filename, d3=False, label=True):
    data = np.loadtxt(filename)
    x = data[:, 0]; y = data[:, 1]; z = data[:, 2]
    if label:
        name = filename.split("/")[-1][:-4]
    else:
        name = None
    if d3:
        ax.plot3D(x, y, z, label=name)
    else:
        ax.plot(x, y, label=name)

def plotSystem(ax, path, d3=False, label=True):
    path = "../data/" + path
    files = glob.glob(path+"/*.dat")
    for filename in files:
        plotPlanet(ax, filename, d3=d3, label=label)


if __name__ == "__main__":
    # sammenligne Euler og Verlet
    figdir = "../figurer/"
    pathEuler = "euler_vs_verlet/euler"
    pathVerlet = "euler_vs_verlet/verlet"

    tlist = [1, 5, 10]
    n = len(tlist)
    fig = plt.figure()
    bigax = emptyax(fig)
    bigax.set_xlabel("x [AU]")
    bigax.set_ylabel("y [AU]")

    labelcounter = 0
    label = True
    for i, t in zip(range(3), tlist):
        if labelcounter > 0:
            label = False
        run_maincpp(1, t, 0.001)
        ax = fig.add_subplot(n, 1, i+1)
        plotSystem(ax, pathEuler, label=label)
        plotSystem(ax, pathVerlet, label=label)
        ax.set_xlim(-1.6, 1.6)
        #plt.axis("equal")
        labelcounter += 1
    fig.legend(fontsize=14)
    plt.savefig(figdir + "eulerVerlet.pdf")


    pathVelo = "escape_velocity"
    pathBeta = "change_beta"

    endtime = 3
    dt = 0.001
    betalist = [2, 2.3, 2.6, 3]
    vscalelist = [1, 1.3, 1.6, 2]
    n = len(betalist)

    vfig = plt.figure()
    vbigax = emptyax(vfig)
    vbigax.set_xlabel("x [AU]")
    vbigax.set_ylabel("y [AU]")

    betafig = plt.figure()
    betabigax = emptyax(betafig)
    betabigax.set_xlabel("x [AU]")
    betabigax.set_ylabel("y [AU]")

    labelcounter = 0
    label = True
    for i, beta, vscale in zip(range(1, 5), betalist, vscalelist):
        if labelcounter > 0:
            label = False
        run_maincpp(2, endtime, dt, vscale, 1)
        ax = vfig.add_subplot(n, 1, i)
        plotSystem(ax, pathVelo, label=label)
        plt.axis("equal")
        run_maincpp(2, endtime, dt, 1, beta)
        ax = betafig.add_subplot(n, 1, i)
        plotSystem(ax, pathBeta, label=label)
        labelcounter += 1
    vfig.legend(fontsize=14)
    betafig.legend(fontsize=14)
    plt.show()

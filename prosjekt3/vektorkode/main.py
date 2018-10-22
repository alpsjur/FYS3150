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


def run_maincpp(n, scenario, endtime, dt, *args):
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

def generatePlots(n, endtime, dt)

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
    if scenario == 2 or scenario == 3:
        velocity_beta = args[0]
        for i in range(n):
            os.system("./main.exe {} {} {} {}".format(scenario, endtime, dt, velocity_beta))
    else:
        for i in range(n):
            os.system("./main.exe {} {} {}".format(scenario, endtime, dt))

def plotPlanet(ax, filename, path=None, d3=False, label=True, centersun=False):
    data = np.loadtxt(filename)
    x = data[:, 0]; y = data[:, 1]; z = data[:, 2]
    if centersun:
        dataS = np.loadtxt(path+"/Sun.dat")
        xS = dataS[:, 0]; yS = dataS[:, 1]; zS = dataS[:, 2]
    else:
        xS = 0; yS=0; zS=0;
    if label:
        name = filename.split("/")[-1][:-4]
    else:
        name = None
    if d3:
        ax.plot3D(x-xS, y-yS, z-yS, label=name)
    else:
        ax.plot(x-xS, y-yS, label=name)

def plotSystem(ax, path, d3=False, label=True, sentersun=False):
    path = "../data/" + path
    files = glob.glob(path+"/*.dat")
    for filename in files:
        plotPlanet(ax, filename, path=path, d3=d3, label=label, centersun = centersun)

def plotAbs(ax, dt, endtime, filename, var, label=None):
    data = np.loadtxt(filename)
    if var == 'pos':
        x = data[:, 0]; y = data[:, 1]; z = data[:, 2]
    if var == 'vel':
        x = data[:, 3]; y = data[:, 4]; z = data[:, 5]
    absVal = np.sqrt(x*x+y*y+z*z)
    t = np.arange(0,endtime,dt)
    ax.plot(t, absVal, label=label)

def plotPosVel(ax, dt, endtime, filename, beta):
    for b in beta:
        run_maincpp(scenario, endtime, dt, b)
        plotAbs(ax[0], dt, endtime, filename, 'pos', label=r"$\beta={}$".format(b))
        plotAbs(ax[1], dt, endtime, filename, 'vel')
        run_maincpp(scenario, endtime, dt, b)
    ax[0].legend(loc='upper center', bbox_to_anchor=(0.5, 1.35),
                  ncol=3,fontsize=12, frameon=False)
    fig.text(0.5, 0.035, 't [yr]',  ha='center',fontsize=14)
    ax[0].set_ylabel('r [AU]')
    ax[1].set_ylabel('v [AU/yr]')

def plotPerihelion(ax, dt, file):
    data = np.loadtxt(file)
    x = data[:, 0]; y = data[:, 1]; i = data[:, 3]
    theta = np.arctan(y/x)
    ax.plot(i*dt,theta)

if __name__ == "__main__":

    sns.set()
    sns.set_style("whitegrid")
    sns.set_palette("husl")
    plt.rc('text', usetex=True)
    figdir = "../figurer/"


    # endrer massen til jupiter
    scenario = 4
    endtime = 20
    dt = 0.001
    masses = [1,10,1000]

    #run_maincpp(scenario, endtime, dt)
    fig, ax = plt.subplots(2, 2)
    sub = [0,2,3]
    for i in range(3):
        plotSystem(ax.flatten()[sub[i]], "sun_earth_jupiter/jupiter_mass_{}".format(masses[i]))
        ax.flatten()[sub[i]].axis('equal')
        i += 1
    ax[1,1].axis([-4, 23, -19, 2])
    ax[1,1].legend(loc='upper center', bbox_to_anchor=(0.5, 2),fontsize=14, frameon=False)
    ax[0,1].axis('off')

    fig.text(0.5, 0.03, 'x [AU]',  ha='center',fontsize=14)
    fig.text(0.02, 0.5, 'y [AU]',  va='center', rotation='vertical',fontsize=14)

    plt.savefig(figdir+"jupiter_mass.pdf")



    #endrer gravitasjonskraften
    scenario = 3
    endtime = [10, 60]
    dt = 0.001
    beta = [2, 2.5, 2.9, 2.99, 2.999, 3]
    filename= "../data/change_beta/Earth.dat"

    fig, ax = plt.subplots(2,1,sharex=True)
    plotPosVel(ax, dt, endtime[0], filename, beta)
    plt.savefig(figdir+"change_beta_10yr.pdf")

    fig, ax = plt.subplots(2,1,sharex=True)
    plotPosVel(ax, dt, endtime[1], filename, beta)
    ax[1].set_ylim(-1,10)
    plt.savefig(figdir+"change_beta_60yr.pdf")

    #massesenter
    scenario = 5
    endtime = 30
    dt = 0.001
    run_maincpp(scenario, endtime, dt)

    fig, ax = plt.subplots(2,2)
    plotSystem(ax[0,1], "sun_earth_jupiter/sun_origo")
    plotSystem(ax[0,0], "sun_earth_jupiter/mass_origo")
    plotSystem(ax[1,1], "sun_earth_jupiter/sun_origo")
    plotSystem(ax[1,0], "sun_earth_jupiter/mass_origo")
    ax[0,0].axis('equal')
    ax[0,1].axis('equal')
    ax[1,0].axis('equal')
    ax[1,1].axis('equal')
    ax[1,1].axis([-0.01,0.07,-0.05,0.003])
    ax[1,0].axis([-0.007,0.007,-0.006,0.006])

    fig.text(0.5, 0.03, 'x [AU]',  ha='center',fontsize=14)
    fig.text(0.02, 0.5, 'y [AU]',  va='center', rotation='vertical', fontsize=14)
    ax[0,0].legend(ncol=3,loc='upper center', bbox_to_anchor=(1, 1.3), fontsize=14,
                    frameon=False)
    #fig.tight_layout()

    #plt.savefig(figdir+"center_of_mass.pdf")



    #hele solsystemet
    scenario = 6
    endtime = 200
    dt = 0.001
    run_maincpp(scenario, endtime, dt)
    sns.set_palette(sns.color_palette("husl", 9))


    fig = plt.figure()
    ax = plt.axes(projection='3d')
    plotSystem(ax, "../data/solarsystem", d3=True)
    ax.legend(loc='upper center',bbox_to_anchor=(0.5, 1.15),
                  ncol=3,fontsize=12, frameon=False)
    ax.set_xlabel('x [AU]')
    ax.set_ylabel('y [AU]')
    ax.set_zlabel('z [AU]')
    #plt.savefig(figdir+"solarsystem3d.pdf")


    fig, ax = plt.subplots()
    plotSystem(ax, "../data/solarsystem")
    ax.legend(loc='upper center',bbox_to_anchor=(0.5, 1.2),
                  ncol=3,fontsize=12, frameon=False)
    ax.axis('equal')
    ax.set_xlabel('x [AU]')
    ax.set_ylabel('y [AU]')
    plt.savefig(figdir+"solarsystem2d.pdf")

    plt.show()

    """
    #The perihelion precession of Mercury
    scenario = 7
    endtime = 1
    dt = 0.0000001

    fig, ax = plt.subplots(1,2)

    plotPerihelion(ax[0], dt, "../data/sun_mercury/classical/MercuryPelihelon.dat")
    plotPerihelion(ax[1], dt, "../data/sun_mercury/relativistic/MercuryPelihelon.dat")
    plt.show()
    """

    """
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
    fig.legend(fontsize=14, frameon=False)
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
    vfig.legend(fontsize=14, frameon=False)
    betafig.legend(fontsize=14, frameon=False)
    plt.show()
    '''

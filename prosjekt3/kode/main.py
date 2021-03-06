import os
import sys
import glob
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np
import seaborn as sns
from scipy.signal import argrelextrema


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

def plotSystem(ax, path, d3=False, label=True, centersun=False):
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
    radToArc = 206264.8062
    data = np.loadtxt(file)
    x = data[:, 0]; y = data[:, 1]; i = data[:, 3]
    theta = np.arctan(y/x)*radToArc
    ax.plot(i*dt,theta)


if __name__ == "__main__":

    sns.set()
    sns.set_style("whitegrid")
    sns.set_palette("husl")
    plt.rc('text', usetex=True)
    figdir = "../figurer/"

    '''
    #unnslipningsfart
    scenario = 2
    endtime = 5
    dt = 0.001
    escape = np.sqrt(8)*np.pi/6.28353
    scale =  [escape*0.9,escape]
    path = "../data/escape_velocity"

    fig, ax = plt.subplots()
    ax.axis("equal")
    for s in scale:
        run_maincpp(scenario, endtime, dt, s)
        plotSystem(ax, path, centersun=True)
    handles = ax.get_legend_handles_labels()
    handle1 = handles[0][1]
    handle2 = handles[0][3]
    ax.legend([handle1,handle2],[r"$0.9 v_e$",r"$v_e$"],fontsize=14)
    ax.set_xlabel("x [AU]",fontsize=14)
    ax.set_ylabel("y [AU]",fontsize=14)
    plt.savefig(figdir+"escape.pdf")




    #timer algoritmene
    scenario = 0
    endtime = [1,10,100,500]
    dt = 0.001
    n = 10
    timeFE = []
    timeVV = []

    for t in endtime:
        if os.path.exists("time{}.log".format(t)):
            os.system("rm time{}.log".format(t))

    for t in endtime:
        for i in range(n):
            os.system("./main.exe {} {} {} >> time{}.log".format(scenario, t, dt, t))


    f = open("../data/timing.dat", "w")
    for i in range(len(endtime)):
        timeFE.append(np.loadtxt("time{}.log".format(endtime[i]))[:,0])
        timeVV.append(np.loadtxt("time{}.log".format(endtime[i]))[:,1])
        relation = np.mean(timeVV[i])/np.mean(timeFE[i])
        error = (np.std(timeFE[i])/np.mean(timeFE[i]) + np.std(timeVV[i])/np.mean(timeVV[i]))*relation
        f.write("{} {:.2E} {:.2E} \n".format(endtime[i],relation,error))


    # endrer massen til jupiter
    scenario = 4
    endtime = 20
    dt = 0.001
    masses = [1,10,1000]

    #run_maincpp(scenario, endtime, dt)
    fig, ax = plt.subplots(2, 2)
    sub = [0,2,3]
    for i in range(3):
        plotSystem(ax.flatten()[sub[i]], "sun_earth_jupiter/jupiter_mass_{}".format(masses[i]),centersun=True)
        ax.flatten()[sub[i]].axis('equal')

    ax[1,1].set_xlim(ax[0,0].get_xlim())
    ax[1,1].set_ylim(ax[0,0].get_ylim())
    ax[1,1].legend(loc='upper center', bbox_to_anchor=(0.5, 2),fontsize=14, frameon=False)
    ax[0,1].axis('off')

    fig.text(0.5, 0.03, 'x [AU]',  ha='center',fontsize=14)
    fig.text(0.02, 0.5, 'y [AU]',  va='center', rotation='vertical',fontsize=14)

    plt.savefig(figdir+"jupiter_mass.pdf")

    plt.show()
    '''
    '''
    #endrer gravitasjonskraften
    scenario = 3
    endtime = [10,60]
    dt = 0.001
    beta = [2, 2.5, 2.9, 2.99, 2.999, 3]
    filename= "../data/change_beta/Earth.dat"

    fig, ax = plt.subplots(2,1,sharex=True)
    plotPosVel(ax, dt, endtime[0], filename, beta)
    #plt.savefig(figdir+"change_beta_10yr.pdf")

    fig, ax = plt.subplots(2,1,sharex=True)
    plotPosVel(ax, dt, endtime[1], filename, beta)
    ax[1].set_ylim(-1,10)
    #plt.savefig(figdir+"change_beta_60yr.pdf")

    '''

    '''
    #massesenter
    scenario = 5
    endtime = 30
    dt = 0.001
    #run_maincpp(scenario, endtime, dt)

    fig = plt.figure()
    ax1 = plt.subplot2grid((2, 2), (0, 0), colspan=1)
    ax2 = plt.subplot2grid((2, 2), (1, 0), colspan=2)
    ax3= plt.subplot2grid((2, 2), (0, 1))

    plotPlanet(ax1, "../data/sun_earth_jupiter/mass_origo/Sun.dat")
    ax1.axis('equal')
    ax1.set_xlabel("x [AU]",fontsize=14)
    ax1.set_ylabel("y [AU]",fontsize=14)

    dataS = np.loadtxt("../data/sun_earth_jupiter/sun_origo/Sun.dat")
    xS = dataS[:, 0]; yS = dataS[:, 1]; zS = dataS[:, 2]

    dataE = np.loadtxt("../data/sun_earth_jupiter/sun_origo/Earth.dat")
    xE = dataE[:, 0]; yE = dataE[:, 1]; zE = dataE[:, 2]

    dataJ = np.loadtxt("../data/sun_earth_jupiter/sun_origo/Jupiter.dat")
    xJ = dataJ[:, 0]; yJ = dataJ[:, 1]; zJ = dataJ[:, 2]

    rE = np.sqrt((xE-xS)**2+ (yE-yS)**2 +(zE-yS)**2)
    rJ = np.sqrt((xJ-xS)**2+ (yJ-yS)**2 +(zJ-yS)**2)

    dataEc = np.loadtxt("../data/sun_earth_jupiter/mass_origo/Earth.dat")
    rEc = np.sqrt(dataEc[:, 0]**2+ dataEc[:, 1]**2 + dataEc[:, 2]**2)
    dataJc = np.loadtxt("../data/sun_earth_jupiter/mass_origo/Jupiter.dat")
    rJc = np.sqrt(dataJc[:, 0]**2+ dataJc[:, 1]**2 + dataJc[:, 2]**2)

    t = np.arange(0,endtime,dt)

    ax2.plot([0,0],[0,0])
    ax2.plot(t,rE-rEc)
    ax2.plot(t,rJ-rJc)
    ax2.set_xlabel("t [yr]",fontsize=14)
    ax2.set_ylabel("correction [AU]",fontsize=14)

    ax3.plot([0,0],label="Sun")
    ax3.plot([0,0],label="Earth")
    ax3.plot([0,0],label="Jupiter")
    ax3.axis([2,3,2,3])
    ax3.legend(loc='upper center', bbox_to_anchor=(0.5, 0.75),fontsize=14, frameon=False)
    ax3.axis('off')

    fig.tight_layout()

    rE = np.sqrt((xE-xS)**2+ (yE-yS)**2 +(zE-yS)**2)
    rJ = np.sqrt((xJ-xS)**2+ (yJ-yS)**2 +(zJ-yS)**2)

    plt.savefig(figdir+"center_of_mass.pdf")

    plt.show()


    #hele solsystemet
    scenario = 6
    endtime = 200
    dt = 0.001
    #run_maincpp(scenario, endtime, dt)
    sns.set_palette(sns.color_palette("husl", 9))


    fig = plt.figure()
    ax = plt.axes(projection='3d')
    plotSystem(ax, "../data/solarsystem", d3=True, centersun=True)
    ax.legend(loc='upper center',bbox_to_anchor=(0.5, 1.15),
                  ncol=3,fontsize=12, frameon=False)
    ax.set_xlabel('x [AU]')
    ax.set_ylabel('y [AU]')
    ax.set_zlabel('z [AU]')
    #plt.savefig(figdir+"solarsystem3d.pdf")


    fig, ax = plt.subplots()
    plotSystem(ax, "../data/solarsystem", centersun=True)
    ax.legend(loc='upper center',bbox_to_anchor=(0.5, 1.2),
                  ncol=3,fontsize=12, frameon=False)
    ax.axis('equal')
    ax.set_xlabel('x [AU]')
    ax.set_ylabel('y [AU]')
    #plt.savefig(figdir+"solarsystem2d.pdf")



    #The perihelion precession of Mercury
    scenario = 7
    endtime = 100
    dt = 0.0000001


    #run_maincpp(scenario, endtime, dt)

    fig, ax = plt.subplots()
    radToArc = 206264.8062
    plotPerihelion(ax, dt, "../data/sun_mercury/classical/MercuryPerihelion.dat")
    plotPerihelion(ax, dt, "../data/sun_mercury/relativistic/MercuryPerihelion.dat")

    ax.set_xlabel("t [yr]",fontsize=14)
    ax.set_ylabel("perihelion precession [arc seconds]",fontsize=14)
    ax.legend(["Classical", "Relativistic"],fontsize=14)

    plt.savefig(figdir+"perihelion.pdf")
    #plt.show()


    # sammenligne Euler og Verlet
    figdir = "../figurer/"
    pathEuler = "euler_vs_verlet/euler"
    pathVerlet = "euler_vs_verlet/verlet"
    tlist = [1, 5, 10]

    fig, ax = plt.subplots(2, 2)
    sub = [0,2,3]
    for i in range(3):
        run_maincpp(1, tlist[i], 0.001)
        plotSystem(ax.flatten()[sub[i]], pathEuler, label=False)
        plotSystem(ax.flatten()[sub[i]], pathVerlet, label=False)
        ax.flatten()[sub[i]].axis('equal')

    ax[0,1].plot([0,0],label=None)
    ax[0,1].plot([0,0],label="Euler")
    ax[0,1].plot([0,0],label=None)
    ax[0,1].plot([0,0],label="Verlet")
    ax[0,1].axis([2,3,2,3])
    ax[0,1].legend(loc='upper center', bbox_to_anchor=(0.5, 0.75),fontsize=14, frameon=False)
    ax[0,1].axis('off')

    fig.text(0.5, 0.03, 'x [AU]',  ha='center',fontsize=14)
    fig.text(0.02, 0.5, 'y [AU]',  va='center', rotation='vertical',fontsize=14)

    plt.savefig(figdir + "eulerVerlet.pdf")
    '''

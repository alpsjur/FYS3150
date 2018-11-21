import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import seaborn as sns


# construcs a simple dataframe from the ising output
def construct_df(filename, header):
    datadir = "../data/"
    df = pd.read_csv(datadir + filename,
                    header=None,
                    delim_whitespace=True,
                    names=header,
                    index_col=0
                    )
    return df

# this function plot all parameters from the ising model
def plot_all(df, dt, header, gridSizes, window=5):
    figdir = "../figurer/"
    # calculating rolling mean of all parameters in the dataframe
    df_rm = df.rolling(window=window,
                        center=True
                        ).mean()
    for parameter in header[1:]:
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        for size in gridSizes:
            ax.scatter(df.index, df[parameter + " {}".format(size)].values)
            df_rm.plot(y=parameter + " {}".format(size), legend=False, label="L={}".format(size), ax=ax)
        ax.set_xlabel(r"Temperature $[\frac{k_B K}{J}]$", fontsize=14)
        ax.set_ylabel(r"{}".format(parameter), fontsize=14)
        fig.legend(loc="upper center", fontsize=14, frameon=False, ncol=4, bbox_to_anchor=[0.5, 1.03])
        plt.tight_layout()
        #plt.savefig(figdir + parameter + dt + ".pdf")

# function to calculate the critical temperature for finite lattices
def calculateTC(df, gridSizes, window=5):
    df_rm = df.rolling(window=window,
                        center=True
                        ).mean()
    tempcrits = np.zeros((2, len(gridSizes)))
    for i in range(len(gridSizes)):
        tempcrits[0, i] = df["Heat capacity {}".format(gridSizes[i])].idxmax(axis=0)
        tempcrits[1, i] = df_rm["Heat capacity {}".format(gridSizes[i])].idxmax(axis=0)
    return tempcrits

# function to estimate the critical temperature in the thermodynamic limit
def estimateTC(tempcrits, gridSizes):
    criticalT = []
    criticalT_rm = []
    for j in range(len(gridSizes)):
        for i in range(j):
            criticalT.append((tempcrits[0, j]*gridSizes[j] - tempcrits[0, i]*gridSizes[i])/(gridSizes[j] - gridSizes[i]))
            criticalT_rm.append((tempcrits[1, j]*gridSizes[j] - tempcrits[1, i]*gridSizes[i])/(gridSizes[j] - gridSizes[i]))
    return np.array([criticalT, criticalT_rm])


# this functions return a empty ax for creating axis labels for multiple subplots
def emptyax(fig):
    ax = fig.add_subplot(1, 1, 1)
    ax.grid(False)
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
    return ax


sns.set()
sns.set_style("whitegrid")
sns.set_palette("husl")
#plt.style.use("seaborn-whitegrid")


gridSizes = [40, 60, 80, 100]
# this section is quite unreadable. header for the ising model
header = ["Temperature", "Energy", "Heat capacity", "Magnetisation",
            "Susceptibility", "Absolute magnetisation"]
# this line creates a list with gridSize behind the parameter names
headertemp = [header[0]] + list(map("".join, zip(header[1:], [" 40" for i in range(len(header)-1)])))
df001 = construct_df("ising_GRID40_1E6_001.dat", headertemp)
df0001 = construct_df("ising_GRID40_1E6_0001.dat", headertemp)

# concatenating all dataframes of all gridsizes
for size in gridSizes[1:]:
    headertemp = [header[0]] + list(map("".join, zip(header[1:], [" {}".format(size) for i in range(len(header)-1)])))
    df001 = pd.concat([df001, construct_df("ising_GRID{}_1E6_001.dat".format(size), headertemp)], axis=1)
    df0001 = pd.concat([df0001, construct_df("ising_GRID{}_1E6_0001.dat".format(size), headertemp)], axis=1)

#print(df001.head())

plot_all(df001, "001", header, gridSizes, window=2)
plot_all(df0001, "0001", header, gridSizes, window=5)
tempcrits001 = calculateTC(df001, gridSizes, window=2)
tempcrits0001 = calculateTC(df0001, gridSizes, window=5)

# haaper framtidige meg ikke ser tilbake paa denne koden
TC001 = estimateTC(tempcrits001, gridSizes)[0, :]
TC001_rm = estimateTC(tempcrits001, gridSizes)[1, :]
TC0001 = estimateTC(tempcrits0001, gridSizes)[0, :]
TC0001_rm = estimateTC(tempcrits0001, gridSizes)[1, :]

print(TC001.mean(), TC001.std())
print(TC001_rm.mean(), TC001_rm.std())
print(TC0001.mean(), TC0001.std())
print(TC0001_rm.mean(), TC0001_rm.std())

figdir = "../figurer/"

# plotting the probability distribution
header = ["energy", "probability", "expval", "STD"]

pdf10 = construct_df("isingPDF_GRID20_1E6_T1.dat", header)
pdf24 = construct_df("isingPDF_GRID20_1E6_T24.dat", header)

#print(pdf10.reset_index().expval, pdf10.reset_index().std)
#print(pdf24.reset_index().expval, pdf24.reset_index().std)

figPDF = plt.figure()

invisible_ax = emptyax(figPDF)          # ax for creating axis labels in the middle
axPDF10 = figPDF.add_subplot(2, 1, 1)
axPDF24 = figPDF.add_subplot(2, 1, 2, sharex=axPDF10)
axPDF10.tick_params(labelbottom=False, bottom=False)

axPDF10.plot(pdf10.index, pdf10.probability, ".")

axPDF24.plot(pdf24.index, pdf24.probability, ".")
axPDF24.set_xlabel("Energy", fontsize=14)
invisible_ax.set_ylabel("Probability", fontsize=14)
figPDF.text(0.75, 0.82, "T = 1.0", fontsize=14)
figPDF.text(0.75, 0.4, "T = 2.4", fontsize=14)
#plt.savefig(figdir + "probabilitydist.pdf")


commstr = "ising_GRID20_MC1E3_25E3_"
mcheader = ["Monte Carlo Cycles", "Energy", "Absolute magnetisation"]

df_T10_ordered = construct_df(commstr + "T10_ordered.dat", mcheader)
df_T10_unordered = construct_df(commstr + "T10_unordered.dat", mcheader)
df_T24_ordered = construct_df(commstr + "T24_ordered.dat", mcheader)
df_T24_unordered = construct_df(commstr + "T24_unordered.dat", mcheader)

fig, axes = plt.subplots(2, 1, sharex=True)
#invisible_ax = emptyax(fig)
for parameter, i in zip(mcheader[1:], range(2)):
    axes[i].plot(df_T10_ordered.index, df_T10_ordered[parameter].values, ".")
    axes[i].plot(df_T10_ordered.index, df_T10_unordered[parameter].values, ".")
    axes[i].plot(df_T24_ordered.index, df_T24_ordered[parameter].values, ".")
    axes[i].plot(df_T24_unordered.index, df_T24_unordered[parameter].values, ".")
    axes[i].set_ylabel(parameter, fontsize=14)
axes[1].set_xlabel("Monte Carlo Cycles", fontsize=14)
#fig.text(0.2, 0.8, "T = 1.0", fontsize=14)
#fig.text(0.2, 0.4, "T = 2.4", fontsize=14)
fig.legend(["ordered T=1.0", "unordered T=1.0", "ordered T=2.4", "unordered T=2.4"],
            loc="upper center",
            frameon=False,
            ncol=2,
            fontsize=14,
            bbox_to_anchor=[0.5, 1.03]
            )
    #ax1.legend(frameon=False)
#plt.savefig(figdir + "MC" + ".pdf")


df_relerror = construct_df("relative_error.dat", mcheader)
df_relerroravg = df_relerror.rolling(window=5,
                    center=True
                    ).mean()
labels = [r"$\epsilon_{E}$", r"$\epsilon_{|M|}$"]

fig, axes = plt.subplots(2, 1, sharex=True)
for i in range(len(mcheader[1:])):
    axes[i].loglog(df_relerror.index, df_relerror[mcheader[i+1]].values, ".")
    axes[i].loglog(df_relerroravg.index, df_relerroravg[mcheader[i+1]].values, linewidth=3, color=sns.color_palette("husl")[2])
    axes[i].set_ylabel(labels[i], fontsize=14)
axes[1].set_xlabel("Monte Carlo Cycles", fontsize=14)
plt.tight_layout()
#plt.savefig(figdir + "relerror.pdf")


df_accept1_unordered = construct_df("acceptance_T1_unordered.dat", ["Monte Carlo cycles", "unordered T=1.0"])
df_accept1_ordered = construct_df("acceptance_T1_ordered.dat", ["Monte Carlo cycles", "ordered T=1.0"])
df_accept24_unordered = construct_df("acceptance_T24_unordered.dat", ["Monte Carlo cycles", "unordered T=2.4"])
df_accept24_ordered = construct_df("acceptance_T24_ordered.dat", ["Monte Carlo cycles", "ordered T=2.4"])

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

df_accept1_unordered.plot(loglog=True, legend=False, ax=ax)
df_accept1_ordered.plot(loglog=True, legend=False, ax=ax)
df_accept24_unordered.plot(loglog=True, legend=False, ax=ax)
df_accept24_ordered.plot(loglog=True, legend=False, ax=ax)

fig.legend(fontsize=14, ncol=2, loc="upper center", frameon=False, bbox_to_anchor=[0.5, 1.03])
ax.set_xlabel("Monte Carlo cycles", fontsize=14)
ax.set_ylabel("Acceptance count", fontsize=14)
#plt.savefig(figdir + "acceptance.pdf")
plt.show()

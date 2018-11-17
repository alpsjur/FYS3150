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
def plot_all(df, header, gridSizes):
    # calculating rolling mean of all parameters in the dataframe
    df_avg = df.rolling(window=4,
                        center=True
                        ).mean()
    print(df_avg.idxmax(axis=0))
    for parameter in header[1:]:
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        for size in gridSizes:
            ax.scatter(df.index, df[parameter + " {}".format(size)].values)
            df_avg.plot(y=parameter + " {}".format(size), legend=False, label="L={}".format(size), ax=ax)
        ax.set_xlabel(r"Temperature $[\frac{K}{J}]$", fontsize=14)
        ax.set_ylabel(r"{}".format(parameter), fontsize=14)
        fig.legend(loc="upper center", fontsize=14, frameon=False, ncol=4)

# this functions return a empty ax for creating axis labels for multiple subplots
def emptyax(fig):
    ax = fig.add_subplot(1, 1, 1)
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
# this line creates a list with _gridSize behind the parameter names
headertemp = [header[0]] + list(map("".join, zip(header[1:], [" 40" for i in range(len(header)-1)])))
df001 = construct_df("ising_GRID40_1E6_001.dat", headertemp)
df0001 = construct_df("ising_GRID40_1E6_0001.dat", headertemp)

# concatenating all dataframes of all gridsizes
for size in gridSizes[1:]:
    headertemp = [header[0]] + list(map("".join, zip(header[1:], [" {}".format(size) for i in range(len(header)-1)])))
    df001 = pd.concat([df001, construct_df("ising_GRID{}_1E6_001.dat".format(size), headertemp)], axis=1)
    df0001 = pd.concat([df0001, construct_df("ising_GRID{}_1E6_0001.dat".format(size), headertemp)], axis=1)

#print(df001.head())

plot_all(df001, header, gridSizes)
plot_all(df0001, header, gridSizes)


# plotting the probability distribution
header = ["energy", "probability"]

pdf10 = construct_df("isingPDF_GRID20_1E6_T1.dat", header)
pdf24 = construct_df("isingPDF_GRID20_1E6_T24.dat", header)

figPDF = plt.figure()

invisible_ax = emptyax(figPDF)          # ax for creating axis labels in the middle
axPDF10 = figPDF.add_subplot(2, 1, 1)
axPDF24 = figPDF.add_subplot(2, 1, 2)

axPDF10.scatter(pdf10.index, pdf10.probability)

axPDF24.plot(pdf24.index, pdf24.probability)
invisible_ax.set_xlabel("Energy", fontsize=14)
invisible_ax.set_ylabel("Probability", fontsize=14)


commstr = "ising_GRID20_MC1E3_25E3_"
mcheader = ["Monte Carlo Cycles", "Energy", "Absolute magnetisation"]

df_T10_ordered = construct_df(commstr + "T10_ordered.dat", mcheader)
df_T10_unordered = construct_df(commstr + "T10_unordered.dat", mcheader)
df_T24_ordered = construct_df(commstr + "T24_ordered.dat", mcheader)
df_T24_unordered = construct_df(commstr + "T24_unordered.dat", mcheader)

figMC = plt.figure()
invisible_ax = emptyax(figMC)
ax1 = figMC.add_subplot(2, 1, 1)
ax2 = figMC.add_subplot(2, 1, 2)

for parameter in mcheader[1:]:
    ax1.plot(df_T10_ordered.index, df_T10_ordered["{}".format(parameter)], label="ordered")
    ax1.plot(df_T10_unordered.index, df_T10_unordered["{}".format(parameter)], label="unordered")
    ax2.plot(df_T24_ordered.index, df_T24_ordered["{}".format(parameter)], label="ordered")
    ax2.plot(df_T24_unordered.index, df_T24_unordered["{}".format(parameter)], label="unordered")
ax1.legend(frameon=False)
ax2.legend(frameon=False)
invisible_ax.set_xlabel("Monte Carlo Cycles")



plt.show()

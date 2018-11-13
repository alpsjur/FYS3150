import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def construct_df(filename):
    datadir = "../data/"
    columns = ["temperature", "energy", "heat capacity", "magnetisation",
                "susceptibility", "something"]
    df = pd.read_csv(datadir + filename,
                    header=None,
                    delim_whitespace=True,
                    names=columns,
                    index_col=0
                    )
    return df


sns.set()
sns.set_style("whitegrid")
sns.set_palette("husl")
plt.rc('text', usetex=True)

grid100 = construct_df("ising_GRID40_1E6_0001.dat")
print(grid100.head())
grid100_avg = grid100.rolling(window=5,
                            #center=True
                            ).mean()
print(grid100_avg.idxmax(axis=0))
#fig = plt.figure()
ax = grid100.reset_index().plot.scatter(x="temperature", y="heat capacity")
grid100_avg.plot(y="heat capacity", ax=ax, legend=False)
ax.set_ylabel("heat capacity [J/kg]", fontsize=14)

plt.show()

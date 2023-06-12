import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import models

def plot_mixed(ind,t,inf,df_seropositives,df_mixed):
    disease = ["Measles","Mumps","Rubella"]

    plt.plot(df_seropositives[0].values,df_seropositives[ind].values,"ko")

    for i in range(n):
        plt.plot(t,models.F(t,*df_mixed.iloc[i].values),label="Outliers: "+str(inf+i))
        plt.title(disease[ind-1])
        plt.legend()

    plt.savefig(disease[ind-1]+".pdf",bbox_inches = "tight")
    plt.show()
    plt.close()

df_seropositives = pd.read_table("output/seropositives.txt",delimiter=" ",header=None,skiprows=1)
df_mixed_measles = pd.read_table("output/solutions_mixed_measles.txt",delimiter=" ",header=None,skiprows=0)
df_mixed_mumps   = pd.read_table("output/solutions_mixed_mumps.txt",delimiter=" ",header=None,skiprows=0)
df_mixed_rubella = pd.read_table("output/solutions_mixed_rubella.txt",delimiter=" ",header=None,skiprows=0)

with open("output/num_mixed_test.txt") as f:
    lines = f.readlines()
    lim = [line.split()[0] for line in lines]

inf = int(lim[0])
sup = int(lim[1])
n = sup - inf + 1
t = np.linspace(0,70,1000)

plot_mixed(1,t,inf,df_seropositives,df_mixed_measles)
plot_mixed(2,t,inf,df_seropositives,df_mixed_mumps)
plot_mixed(3,t,inf,df_seropositives,df_mixed_rubella)



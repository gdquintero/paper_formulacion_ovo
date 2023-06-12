import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import models

# 'ro',mfc='none',ms=10)

def plot_solutions(ind,df_seropositives,df_sol_ovo1,df_sol_ovo2,df_sol_ls):
    t = np.linspace(0,70,1000)
    disease = ["Measles","Mumps","Rubella"]
    plt.plot(df_seropositives[0].values,df_seropositives[ind].values,"ko")
    # plt.plot(t,models.F(t,*solutions_farrington[i-1]),label="Farrington")
    plt.plot(t,models.F(t,*df_sol_ls.iloc[ind-1].values),label="Least Squares")
    plt.plot(t,models.F(t,*df_sol_ovo1.iloc[ind-1].values),label="OVO with initial point")
    plt.plot(t,models.F(t,*df_sol_ovo2.iloc[ind-1].values),label="OVO with midpoint")

    plt.legend()
    plt.title(disease[ind-1])
    plt.savefig(disease[ind-1]+".pdf",bbox_inches = "tight")
    # plt.show()
    plt.close()

df_seropositives = pd.read_table("output/seropositives.txt",delimiter=" ",header=None,skiprows=1)
df_solutions_t1 = pd.read_table("output/solutions_ovo_wins_t1.txt",delimiter=" ",header=None,skiprows=0)
df_solutions_t2 = pd.read_table("output/solutions_ovo_wins_t2.txt",delimiter=" ",header=None,skiprows=0)
df_solutions_ls  = pd.read_table("output/solutions_ls.txt",delimiter=" ",header=None,skiprows=0)


# Plotamos las soluciones 1:Measles, 2:Mumps, 3:Rubella
plot_solutions(1,df_seropositives,df_solutions_t1,df_solutions_t2,df_solutions_ls)
plot_solutions(2,df_seropositives,df_solutions_t1,df_solutions_t2,df_solutions_ls)
plot_solutions(3,df_seropositives,df_solutions_t1,df_solutions_t2,df_solutions_ls)


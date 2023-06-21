import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

model = lambda t,x1,x2,x3,x4 : x1 + x2*t + x3*(t**2) + x4*(t**3)

xstar = np.array([0,2,-3,1])

t = np.linspace(-1.1,3.6,1000)

# df_solutions =  pd.read_table("output/solutions.txt",delimiter=" ",header=None,skiprows=0)
df_data = pd.read_table("output/data.txt",delimiter=" ",header=None,skiprows=1)

with open("output/solutions.txt") as f:
    lines = f.readlines()
    solutions = [line.split()[:] for line in lines]

with open("output/num_mixed_test.txt") as f:
    lines = f.readlines()
    lim = [line.split()[0] for line in lines]


inf = int(lim[0])
sup = int(lim[1])
n = sup - inf + 1

plt.plot(df_data[0].values,df_data[1].values,"ko")

for i in range(n):
    plt.plot(t,model(t,float(solutions[i][0]),float(solutions[i][1]),float(solutions[i][2]),\
                     float(solutions[i][3])),label="Outliers: "+str(inf+i))
    plt.legend()
    

plt.show()
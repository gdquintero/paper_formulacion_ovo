import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random

model = lambda t,x1,x2,x3,x4 : x1 + x2*t + x3*(t**2) + x4*(t**3)

xstar = np.array([0,2,-3,1])

m = 46

t = np.linspace(-1,3.5,m)

y = 10 * np.ones(m)

for i in range(m):
    if i <= 5 or 16 <= i:
        y[i] = model(t[i],*xstar) + random.uniform(-0.01,0.01)

with open("output/data.txt","w") as f:
    f.write("%i\n" % m)
    for i in range(m):
        f.write("%f %f\n" % (t[i],y[i]))

plt.plot(t,y,"o")
plt.show()
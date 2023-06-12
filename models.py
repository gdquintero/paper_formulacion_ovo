import numpy as np

l = lambda t, a, b, c: (a * t - c) * np.exp(-b * t) + c

F = lambda t, a, b, c: 1.0 - np.exp((a/b) * t * np.exp(-b * t) + (1.0/b) * ((a/b) - c) * \
    (np.exp(-b * t) - 1.0) - c * t)
#!/anaconda3/bin/python

import numpy as np
import sys
sys.path.append("./wrapper")
sys.path.append("./mesh")
from wrapper_flux import c_roe_2d, c_flux_function_2d, c_roe_2d


uL = np.array([1.0, 2.0, 3.0, 17.0])
uR = np.array([2.0, 3.0, 4.0, 50.0])
n = np.array([1, 2] / np.sqrt(5))
gamma = 1.4

print("=================================")
print(c_roe_2d(uL, uL, n, gamma))
print(c_flux_function_2d(uL, gamma))
print(np.matmul(c_flux_function_2d(uL, gamma), n))

print(c_roe_2d(uL, uR, n, gamma)[0])
print(c_roe_2d(uR, uL, -1.0 * n, gamma)[0])

print("=================================")
uL = np.array([1.0, 10.0, 0.0, 20.0])
uR = np.array([1.0, 1.0, 0.0, 20.0])
n = np.array([1, 0])
gamma = 1.4
print(np.reshape(c_roe_2d(uL, uR, n, gamma)[0], (1, 4)))
print(np.matmul(c_flux_function_2d(uL, gamma), n))
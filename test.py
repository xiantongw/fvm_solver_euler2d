#!/anaconda3/bin/python

import numpy as np
import time
import sys
sys.path.append("./wrapper")
sys.path.append("./mesh")

from wrapper_flux import c_roe_2d, c_flux_function_2d, c_roe_2d
from wrapper_euler_steady_solver import c_euler_solver_main
from TriMesh import TriMesh


uL = np.array([2.0, 2.0, 2.0, 10.0])
uR = np.array([2.0, 1.0, 1.0, 12.0])
n = np.array([np.sqrt(2)/2, np.sqrt(2)/2])
gamma = 1.4

#print(c_roe_2d(uL, uR, n, gamma))
# print(c_roe_2d(uR, uL, -1.0 * n, gamma))

# print(c_roe_2d(uL, uL, n, gamma))
# print(np.matmul(c_flux_function_2d(uL, gamma), n))

mesh = TriMesh("./mesh/bump3.gri")
mesh.info()
mesh.verification()
print(mesh.Bname)

print("****************************************")
print("           Simulation Start             ")
print("****************************************")

param = {
    "cfl": 0.85,
    "mach_inf": 0.5,
    "attack_angle": 0,
    "gamma": 1.4,
    "p_inf": 1.0,
    "bound0":"Inviscid_Wall",
    "bound1":"Subsonic_Outflow",
    "bound2":"Inviscid_Wall",
    "bound3":"Inflow",
    "MAXITER": 100000,
    "eps": 1e-8,
    "nstage":2}

c_euler_solver_main(mesh, param)

    

#!/anaconda3/bin/python

import numpy as np
import time
import sys
sys.path.append("./wrapper")
sys.path.append("./mesh")

from wrapper_flux import c_roe_2d, c_flux_function_2d, c_roe_2d
from wrapper_euler_steady_solver import c_euler_solver_main
from TriMesh import TriMesh
from plot_state import plot_state

uL = np.array([1.0, 2.0, 3.0, 17.0])
uR = np.array([2.0, 3.0, 4.0, 50.0])
n = np.array([1, 2])
gamma = 1.4

#print(c_roe_2d(uL, uR, n, gamma))
# print(c_roe_2d(uR, uL, -1.0 * n, gamma))

# print(c_roe_2d(uL, uL, n, gamma))
# print(np.matmul(c_flux_function_2d(uL, gamma), n))

mesh = TriMesh("./mesh/bump3.gri")


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
 "eps": 1e-7,
 "nstage":2
}

    
state_vectors = c_euler_solver_main(mesh, param)

data = state_vectors[:, 0]

fig, ax = plot_state(mesh, data, show_mesh=True,
                     x_min=-1.5, x_max=1.5, y_min=0, y_max=0.8)

fig.savefig("test.eps")




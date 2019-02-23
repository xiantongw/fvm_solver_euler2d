#!/anaconda3/bin/python

import numpy as np
import time
import sys
sys.path.append("./wrapper")
sys.path.append("./mesh")

from wrapper_flux import c_roe_2d, c_flux_function_2d, c_roe_2d
from wrapper_euler_steady_solver import c_euler_solver_main
from TriMesh import TriMesh


mesh = TriMesh("./mesh/bump3.gri")


param = {
 "cfl": 0.85,
 "mach_inf": 0.50,
 "attack_angle": 0,
 "gamma": 1.4,
 "p_inf": 1.0,
 "bound0": "Inviscid_Wall",
 "bound1": "Subsonic_Outflow",
 "bound2": "Inviscid_Wall",
 "bound3": "Inflow",
 "MAXITER": 1000000,
 "eps": 1e-7,
 "norder": 1,
 }
    
state_vectors = c_euler_solver_main(mesh, param)
np.save("states_bump3_1st.npy", state_vectors)




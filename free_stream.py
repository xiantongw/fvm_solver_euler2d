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


mesh = TriMesh("./mesh/bump4.gri")


param = {
 "cfl": 0.85,
 "mach_inf": 0.50,
 "attack_angle": 0,
 "gamma": 1.4,
 "p_inf": 1.0,
 "bound0":"free_stream",
 "bound1":"free_stream",
 "bound2":"free_stream",
 "bound3":"free_stream",
 "MAXITER": 1000,
 "eps": 1e-20,
 "norder":1
}

    
state_vectors = c_euler_solver_main(mesh, param)





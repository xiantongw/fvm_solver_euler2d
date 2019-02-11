#!/anaconda3/bin/python

import numpy as np
import time
import sys
sys.path.append("./wrapper")
sys.path.append("./mesh")

from wrapper_flux import c_roe_2d, c_flux_function_2d, c_roe_2d
from wrapper_euler_steady_solver import c_euler_solver_main
from TriMesh import TriMesh

def plot_state(mesh, data, show_mesh=True,
               x_min=-1.5, x_max=1.5, y_min=0, y_max=0.8,
               cmap='jet', title='Param_title',
               xlabel='Param_xlabel', ylabel='Param_ylabel',
               cbar_title='Param_cbar_title'):
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    fig, ax = plt.subplots()
    ax.set_aspect('equal', 'box')
    ax.set_xlim([x_min, x_max])
    ax.set_ylim([y_min, y_max])
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if show_mesh:
        cntr = ax.tripcolor(mesh.V[:, 0], mesh.V[:, 1], mesh.E - 1,
                          facecolors=data, cmap=cmap, edgecolors='k', linewidth=0.01)
    else:
        cntr = ax.tripcolor(mesh.V[:, 0], mesh.V[:, 1], mesh.E - 1,
                          facecolors=data, cmap=cmap)
    # create an axes on the right side of ax. The width of cax will be 2%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.08)
    cbar = plt.colorbar(cntr, cax=cax)
    cbar.set_label(cbar_title)
    return fig, ax



uL = np.array([1.0, 2.0, 3.0, 17.0])
uR = np.array([2.0, 3.0, 4.0, 50.0])
n = np.array([1, 2])
gamma = 1.4

print(c_roe_2d(uL, uR, n, gamma))
# print(c_roe_2d(uR, uL, -1.0 * n, gamma))

# print(c_roe_2d(uL, uL, n, gamma))
# print(np.matmul(c_flux_function_2d(uL, gamma), n))

mesh = TriMesh("./mesh/bump0.gri")
mesh.info()
mesh.verification()

# print(mesh.E)
# print(mesh.V)
# print(mesh.I2E)
# print(mesh.B2E)
# print(mesh.In)
# print(mesh.Bn)

print("****************************************")
print("           Simulation Start             ")
print("****************************************")

#import matplotlib
#matplotlib.use('Agg')

for niter in [5000]:

    do_freestream = False
    
    if not do_freestream:
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
         "MAXITER": 1000000,
         "eps": 1e-6,
         "nstage":2}
        
    else:
        param = {
        "cfl": 0.4,
        "mach_inf": 0.5,
        "attack_angle": 0.0,
        "gamma": 1.4,
        "p_inf": 1.0,
        "bound0":"free_stream",
        "bound1":"free_stream",
        "bound2":"free_stream",
        "bound3":"free_stream",
        "MAXITER": niter,
        "eps": 1e-20,
        "nstage":2}
    
    state_vectors = c_euler_solver_main(mesh, param)
    
    #print(residual)

    # print(np.max(gradu), np.min(gradu))
    # print(np.max(residual), np.min(residual))
    
    data = state_vectors[:, 2]

    fig, ax = plot_state(mesh, data, show_mesh=False,
                         x_min=-1.5, x_max=1.5, y_min=0, y_max=0.8)
   
    fig.savefig("{:03d}".format(niter)+".eps")





#test_states = np.zeros([np.shape(mesh.E)[0], 4])

# for i in range(np.shape(mesh.E)[0]):
#         test_states[i][0] = 1.0
#         test_states[i][1] = param['mach_inf'] * np.cos(param['attack_angle'])
#         test_states[i][2] = param['mach_inf'] * np.sin(param['attack_angle'])
#         test_states[i][3] = 1 / (param['gamma'] * (param['gamma'] - 1)) + 0.5 * param['mach_inf'] * param['mach_inf']

#gradu = c_calculate_gradient(mesh, param, test_states)

#print(gradu)


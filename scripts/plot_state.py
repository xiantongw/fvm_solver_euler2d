import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from TriMesh import TriMesh
import numpy as np

def plot_state(mesh, data, show_mesh=False,
               x_min=-1.5, x_max=1.5, y_min=0, y_max=0.8,
               cmap='viridis', title='Param_title',
               xlabel=r'$x$', ylabel=r'$y$',
               cbar_title='Param_cbar_title'):

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

mesh_name = 'bump4'
name_var = 'rho'
mesh = TriMesh("./mesh/{:s}.gri".format(mesh_name))
u = np.loadtxt("./build/states_m0.85_{:s}.dat".format(mesh_name))
gamma = 1.40; p_inf = 1.0; M_inf = 3.0
p = (gamma - 1.0) * (u[:, 3] - 0.5 * (u[:, 1]*u[:, 1] +
            u[:, 2]*u[:, 2]) / u[:, 0])
c = np.sqrt(gamma * p / u[:, 0])
v = np.sqrt((np.power(u[:, 1], 2) + np.power(u[:, 2], 2)) 
                / np.power(u[:, 0], 2))
mach = v/c
states = {'mach':mach, 'p':p, 'ux':u[:, 1] / u[:, 0], 'uy':u[:, 2] / u[:, 0], 'rho':u[:, 0]}
string_title = {'mach':r'$M$', 'p':r'$p$', 'ux':r'$u_x$', 'uy':r'$u_y$', 'rho':r'$\rho$'}


data = states[name_var]
fig, ax = plot_state(mesh, data, cmap='jet', title=string_title[name_var], cbar_title=string_title[name_var])
fig.savefig("{:s}_{:s}_m0.85.eps".format(name_var, mesh_name), bbox_inches='tight')
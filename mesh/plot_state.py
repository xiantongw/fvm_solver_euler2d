import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
    
def plot_state(mesh, data, show_mesh=True,
               x_min=-1.5, x_max=1.5, y_min=0, y_max=0.8,
               cmap='jet', title='Param_title',
               xlabel='Param_xlabel', ylabel='Param_ylabel',
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
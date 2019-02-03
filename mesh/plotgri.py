# -*- coding: utf-8 -*-

from TriMesh import TriMesh
import matplotlib.pyplot as plt

def plotgri(grifile):
    mesh = TriMesh(grifile)
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    ax.triplot(mesh.V[:, 0], mesh.V[:, 1], mesh.E - 1, lw=0.1)
    return fig, ax


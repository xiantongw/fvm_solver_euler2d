#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 13 17:32:45 2019

@author: xtwang
"""
import gmsh
import numpy as np
import meshio
from bump_func import bump_func

def gen_mesh(filename, lc = 0.5, n_bump_points = 21):
    # Set up boundary points and geometry
    x = np.linspace(-0.5,0.5,n_bump_points)
    y = bump_func(x)
    x = np.append(x, [1.5, 1.5, -1.5, -1.5])
    y = np.append(y, [0.0, 0.8, 0.8, 0])
    # Set up GMSH
    print('==================== GMSH START ====================')
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("Bump")
    for iPoint in range(len(x)):
        gmsh.model.geo.addPoint(x[iPoint], y[iPoint], 0.0, lc, iPoint+1)
    for iPoint in range(len(x) - 1):
        gmsh.model.geo.addLine(iPoint+1, iPoint+2)
    gmsh.model.geo.addLine(len(x), 1)
    gmsh.model.geo.addCurveLoop(np.arange(len(x)) + 1, 1)
    gmsh.model.geo.addPlaneSurface([1], 1)
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)
    gmsh.write("bump.msh")
    gmsh.finalize()
    print('==================== GMSH END ======================')
    # Finished Mesh generation
    
    # Prosess the mesh information
    print('Processing gri ')
    mesh_msh = meshio.read('bump.msh')
    E = mesh_msh.cells['triangle']
    V = mesh_msh.points[:, 0:2]
    B_total = mesh_msh.cells['line']
    
    Bname = ['Bottom', 'Right', 'Top', 'Left']
    # temp arrays for 4 boundaries
    B1 = np.empty((0, 2), int)
    B2 = np.empty((0, 2), int)
    B3 = np.empty((0, 2), int)
    B4 = np.empty((0, 2), int)
    # Examine each boundary edge
    for i in range(len(B_total)):
        p1 = V[B_total[i, 0]]
        p2 = V[B_total[i, 1]]
        if abs(p1[0] - 1.5) < 1e-6 and abs(p2[0] - 1.5) < 1e-6:
            B2 = np.vstack([B2, B_total[i, :]])
        elif abs(p1[1] - 0.8) < 1e-6 and abs(p2[1] - 0.8) < 1e-6:
            B3 = np.vstack([B3, B_total[i, :]])
        elif abs(p1[0] + 1.5) < 1e-6 and abs(p2[0] + 1.5) < 1e-6:
            B4 = np.vstack([B4, B_total[i, :]])
        else:
            B1 = np.vstack([B1, B_total[i, :]])
    
    B = [B1, B2, B3, B4]
    print("Processing gri...Done")
    print("=================================================")
    print("Mesh Info:")
    print("Number of Elements: {:d}".format(len(E)))
    print("Number of Nodes: {:d}".format(len(V)))
    print("Writing to File")
    f = open(filename,"w")
    f.write("{:d} {:d} {:d}\n".format(np.shape(V)[0], np.shape(E)[0], np.shape(V)[1]))
    for i in range(np.shape(V)[0]):
        f.write("{:.8E} {:.8E}\n".format(V[i, 0], V[i, 1]))
    f.write("{:d}\n".format(len(Bname)))
    for iB in range(len(Bname)):
        # write the bounday info (1 line)
        f.write("{:d} {:d} {:s}\n".format(len(B[iB]), 2, Bname[iB]))
        for j in range(len(B[iB])):
            f.write("{:d} {:d}\n".format(B[iB][j, 0] + 1, B[iB][j, 1] + 1))
    f.write("{:d} {:d} {:s}\n".format(len(E), 1, "TriLagrange"))
    for i in range(len(E)):
        f.write("{:d} {:d} {:d}\n".format(E[i, 0] + 1, E[i, 1] + 1, E[i, 2] + 1))
    f.close()
    print("Writing to File...Done Filename: "+filename)
    
    
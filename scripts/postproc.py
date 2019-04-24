#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 17 19:14:44 2019

@author: xtwang
"""


import sys
sys.path.append("./wrapper")
sys.path.append("./mesh")
from TriMesh import TriMesh
import numpy as np
import matplotlib.pyplot as plt

def postproc(mesh, param, u):
    gamma = param['gamma']
    p_inf = param['p_inf']
    M_inf = param['mach_inf']
    h = 0.0625
    norder = param['norder']
    # calculate the pressure of each cell
    p = (gamma - 1.0) * (u[:, 3] - 0.5 * (u[:, 1]*u[:, 1] +
            u[:, 2]*u[:, 2]) / u[:, 0])
    s = p / np.power(u[:, 0], gamma)
    
    # calculate the stagnation entropy
    T_t = 1.0 + 0.5 * (gamma - 1) * M_inf ** 2
    p_t = np.power(T_t, gamma / (gamma - 1))
    R = 1.0
    rho_t = p_t / (R * T_t)
    s_t = p_t / (rho_t ** gamma)
    Es = np.sqrt(np.dot((s/s_t - 1.0)**2, mesh.Area) 
                        / np.sum(mesh.Area))[0]
    ind_bottom = np.where(mesh.B2E[:,2] == 1)[0]
    ind_bottom_cells = mesh.B2E[ind_bottom][:, 0] - 1
    
    intg1 = 0.0
    intg2 = 0.0
    cp = [] # pressure coefficient
    mach = [] # mach number
    
    for i in range(len(ind_bottom)):
        # calculate the dl of the bottom cell
        loc_ind = np.array(list({1, 2, 3} - {mesh.B2E[i, 1]})) - 1
        indA, indB = mesh.E[ind_bottom_cells[i], loc_ind] - 1
        dl = np.linalg.norm(mesh.V[indA, :] - mesh.V[indB, :])
        mid_point = 0.5 * (mesh.V[indA, :] + mesh.V[indB, :])
        if norder == 1:
            # first order scheme, using pressure at the cell center
            ul = u[ind_bottom_cells[i], :]
            
        pl = (gamma - 1.0) * (ul[3] - 0.5 * (ul[1]*ul[1] +
                     ul[2]*ul[2]) / ul[0])
        c = np.sqrt(gamma * pl / ul[0])
        m = np.sqrt((ul[1] / ul[0])**2 + (ul[2] / ul[0])**2) / c
        
        intg1 += (pl - p_inf) * mesh.Bn[ind_bottom[i], 1] * dl
        intg2 += (pl - p_inf) * mesh.Bn[ind_bottom[i], 0] * dl
        cp.append([mid_point[0], pl - p_inf])
        mach.append([mid_point[0], m])
        
    cl = intg1 / (gamma * p_inf * (M_inf**2) * h / 2)
    cd = intg2 / (gamma * p_inf * (M_inf**2) * h / 2)
    cp = np.array(cp) / (gamma * p_inf * (M_inf**2) / 2)
    mach = np.array(mach)
    
    cp = cp[cp[:, 0].argsort()]
    mach = mach[mach[:, 0].argsort()]
    
    return cl, cd, cp, Es, mach

param = {
 "cfl": 0.85,
 "mach_inf": 0.85,
 "attack_angle": 0,
 "gamma": 1.4,
 "p_inf": 1.0,
 "bound0":"Inviscid_Wall",
 "bound1":"Subsonic_Outflow",
 "bound2":"Inviscid_Wall",
 "bound3":"Inflow",
 "MAXITER": 100000,
 "eps": 1e-7,
 "norder":1
}


mesh_name = 'bump0'
name_var = 'mach'
mesh = TriMesh("./mesh/{:s}.gri".format(mesh_name))
u = np.loadtxt("./build/states_m0.85_{:s}.dat".format(mesh_name))
cl, cd, cp, Es, mach = postproc(mesh, param, u)
plt.plot(cp[:,0], -cp[:,1],'.-')
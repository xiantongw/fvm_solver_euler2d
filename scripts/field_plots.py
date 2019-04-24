#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 17:08:02 2019

@author: xtwang
"""
import numpy as np
import sys
sys.path.append("./wrapper/")
sys.path.append("./mesh/")
from TriMesh import TriMesh
from plot_state import plot_state

def cal_mach_p(mesh, param, u):
    gamma = param['gamma']
    p_inf = param['p_inf']
    M_inf = param['mach_inf']
    p = (gamma - 1.0) * (u[:, 3] - 0.5 * (u[:, 1]*u[:, 1] +
            u[:, 2]*u[:, 2]) / u[:, 0])
    cp = (p - p_inf) / (0.5 * gamma * p_inf * M_inf * M_inf)
    c = np.sqrt(gamma * p / u[:, 0])
    v = np.sqrt((np.power(u[:, 1], 2) + np.power(u[:, 2], 2)) 
                    / np.power(u[:, 0], 2))
    mach = v/c
    return cp, mach




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
 "norder":2
}
u_1 = np.load("states_bump2_1st.npy")
u_2 = np.load("states_bump2_2nd.npy")
mesh = TriMesh('./mesh/bump2.gri')


cp_1, mach_1 = cal_mach_p(mesh, param, u_1)
fig, ax = plot_state(mesh, cp_1, cmap='jet', title='Pressure Coefficient - 1st order scheme',
           xlabel='x', ylabel='y',
           cbar_title=r'$c_p$')
fig.savefig('cp_1st.eps', bbox_inches='tight')

fig, ax = plot_state(mesh, mach_1, cmap='jet', title='Mach Number - 1st order scheme',
           xlabel='x', ylabel='y',
           cbar_title='Mach Number')
fig.savefig('mach_1st.eps', bbox_inches='tight')

cp_2, mach_2 = cal_mach_p(mesh, param, u_2)
fig, ax = plot_state(mesh, cp_2, cmap='jet', title='Pressure Coefficient - 2nd order scheme',
           xlabel='x', ylabel='y',
           cbar_title=r'$c_p$')
fig.savefig('cp_2nd.eps', bbox_inches='tight')
fig, ax = plot_state(mesh, mach_2, cmap='jet', title='Mach Number - 2nd order scheme',
           xlabel='x', ylabel='y',
           cbar_title='Mach Number')
fig.savefig('mach_2nd.eps', bbox_inches='tight')
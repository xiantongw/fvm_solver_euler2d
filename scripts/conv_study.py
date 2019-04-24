#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 16:24:22 2019

@author: xtwang
"""
import numpy as np
import matplotlib.pyplot as plt
from postproc import postproc
import sys
sys.path.append("./wrapper/")
sys.path.append("./mesh/")
from TriMesh import TriMesh



cl_exact = 1.537095
cd_exact = 2.94278e-6
Es_exact = 0.0

meshes = ['bump0', 'bump1', 'bump2', 'bump3', 'bump4']
#meshes = ['bump0', 'bump1', 'bump2']

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

err_cl_1 = []; err_cd_1 = []; err_Es_1 = []; sqrt_dof=[]
err_cl_2 = []; err_cd_2 = []; err_Es_2 = []
# convergence study on first order scheme
for mesh_name in meshes:
    u = np.load("states_{:s}_1st.npy".format(mesh_name))
    mesh = TriMesh("./mesh/{:s}.gri".format(mesh_name))
    cl, cd, cp, Es, mach = postproc(mesh, param, u)
    dof = len(mesh.E)
    sqrt_dof.append(np.sqrt(dof))
    err_cl_1.append(np.abs(cl - cl_exact))
    err_cd_1.append(np.abs(cd - cd_exact))
    err_Es_1.append(np.abs(Es - Es_exact))

for mesh_name in meshes:
    u = np.load("states_{:s}_2nd.npy".format(mesh_name))
    mesh = TriMesh("./mesh/{:s}.gri".format(mesh_name))
    cl, cd, cp, Es, mach = postproc(mesh, param, u)
    err_cl_2.append(np.abs(cl - cl_exact))
    err_cd_2.append(np.abs(cd - cd_exact))
    err_Es_2.append(np.abs(Es - Es_exact))

fig = plt.figure(figsize=(20,4))

ax1 = fig.add_subplot(1, 3 ,1)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel(r'$\sqrt{dof}$')
ax1.set_ylabel('err')
ax1.set_title('Error of Lift Coefficient')
r1 = np.polyfit(np.log(sqrt_dof), np.log(err_cl_1), 1)[0]
r2 = np.polyfit(np.log(sqrt_dof), np.log(err_cl_2), 1)[0]
ax1.plot(sqrt_dof, err_cl_1, marker='s',
             label='1st order, rate={:.5f}'.format(r1))
ax1.plot(sqrt_dof, err_cl_2, marker='s',
             label='2nd order, rate={:.5f}'.format(r2))
plt.legend()

ax2 = fig.add_subplot(1, 3 ,2)
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel(r'$\sqrt{dof}$')
ax2.set_ylabel('err')
ax2.set_title('Error of Drag Coefficient')
r1 = np.polyfit(np.log(sqrt_dof), np.log(err_cd_1), 1)[0]
r2 = np.polyfit(np.log(sqrt_dof), np.log(err_cd_2), 1)[0]
ax2.plot(sqrt_dof, err_cd_1, marker='s',
             label='1st order, rate={:.5f}'.format(r1))
ax2.plot(sqrt_dof, err_cd_2, marker='s',
             label='2nd order, rate={:.5f}'.format(r2))
plt.legend()

ax3 = fig.add_subplot(1, 3 ,3)
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.set_xlabel(r'$\sqrt{dof}$')
ax3.set_ylabel('err')
ax3.set_title('Error of Entropy Error')
r1 = np.polyfit(np.log(sqrt_dof), np.log(err_Es_1), 1)[0]
r2 = np.polyfit(np.log(sqrt_dof), np.log(err_Es_2), 1)[0]
ax3.plot(sqrt_dof, err_Es_1, marker='s',
             label='1st order, rate={:.5f}'.format(r1))
ax3.plot(sqrt_dof, err_Es_2, marker='s',
             label='2nd order, rate={:.5f}'.format(r2))
plt.legend()

plt.savefig('err_conv.eps', bbox_inches='tight')

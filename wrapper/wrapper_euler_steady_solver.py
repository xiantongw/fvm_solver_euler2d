from ctypes import *
import numpy as np

def c_euler_solver_main(mesh, param):
    solver = np.ctypeslib.load_library("libs/euler_steady_solver.so", ".")

    class CfdParam(Structure):
        _fields_ = [("cfl", c_double),
                    ("mach_inf", c_double),
                    ("attack_angle", c_double),
                    ("gamma", c_double),
                    ("p_inf", c_double)]

    class MeshParam(Structure):
        _fields_ = [("nelem", c_int),
                    ("nnode", c_int),
                    ("niedge", c_int),
                    ("nbedge", c_int),
                    ("nbgroup", c_int)]

    class BoundaryParam(Structure):
        _fields_ = [("bound0", c_char_p),
                    ("bound1", c_char_p),
                    ("bound2", c_char_p),
                    ("bound3", c_char_p)]

    cfd_param = CfdParam(param['cfl'], param['mach_inf'], 
                        param['attack_angle'], param['gamma'], param["p_inf"])

    mesh_param = MeshParam(np.shape(mesh.E)[0], np.shape(mesh.V)[0], np.shape(mesh.I2E)[0],
                            np.shape(mesh.B2E)[0], np.shape(mesh.Bname)[0])
    
    str_bound0 = param["bound0"].encode('utf-8')
    str_bound1 = param["bound1"].encode('utf-8')
    str_bound2 = param["bound2"].encode('utf-8')
    str_bound3 = param["bound3"].encode('utf-8')
    boudary_param = BoundaryParam(str_bound0, str_bound1, str_bound2, str_bound3)

    # Passing arrays from numpy array to C
    c_E = ((c_int * 3) * np.shape(mesh.E)[0])()
    c_V = ((c_double * 2) * np.shape(mesh.V)[0])()
    c_I2E = ((c_int * 4) * np.shape(mesh.I2E)[0])()
    c_B2E = ((c_int * 3) * np.shape(mesh.B2E)[0])()
    c_In = ((c_double * 2) * np.shape(mesh.In)[0])()
    c_Bn = ((c_double * 2) * np.shape(mesh.Bn)[0])()
    c_Area = (c_double * np.shape(mesh.Area)[0])()
    c_state_vectors = ((c_double * 4) * np.shape(mesh.E)[0])()

    # Set up the free strem initial states
    for i in range(np.shape(mesh.E)[0]):
        c_state_vectors[i][0] = 1.0
        c_state_vectors[i][1] = param['mach_inf'] * np.cos(param['attack_angle'])
        c_state_vectors[i][2] = param['mach_inf'] * np.sin(param['attack_angle'])
        c_state_vectors[i][3] = 1 / (param['gamma'] * (param['gamma'] - 1)) + 0.5 * param['mach_inf'] * param['mach_inf']


    for i in range(np.shape(mesh.E)[0]):
        for j in range(3):
            c_E[i][j] = mesh.E[i][j]

    for i in range(np.shape(mesh.V)[0]):
        for j in range(2):
            c_V[i][j] = mesh.V[i][j]

    for i in range(np.shape(mesh.I2E)[0]):
        for j in range(4):
            c_I2E[i][j] = mesh.I2E[i][j]

    for i in range(np.shape(mesh.B2E)[0]):
        for j in range(3):
            c_B2E[i][j] = mesh.B2E[i][j]
    
    for i in range(np.shape(mesh.In)[0]):
        for j in range(2):
            c_In[i][j] = mesh.In[i][j]

    for i in range(np.shape(mesh.Bn)[0]):
        for j in range(2):
            c_Bn[i][j] = mesh.Bn[i][j]

    for i in range(np.shape(mesh.Area)[0]):
        c_Area[i] = mesh.Area[i]
    
    stage = param["nstage"]
    if stage == 1:
        solver.euler_solver_first_order(byref(cfd_param), byref(mesh_param), byref(boudary_param),
                            c_E, c_V, c_I2E, c_B2E, c_In, c_Bn, c_Area,
                            c_state_vectors,
                            c_int(param["MAXITER"]), c_double(param["eps"]))
    if stage == 2:
        solver.euler_solver_second_order(byref(cfd_param), byref(mesh_param), byref(boudary_param),
                            c_E, c_V, c_I2E, c_B2E, c_In, c_Bn, c_Area,
                            c_state_vectors,
                            c_int(param["MAXITER"]), c_double(param["eps"]))


    

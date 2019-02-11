from ctypes import *
import numpy as np

def c_calculate_gradient(mesh, param, state_vectors):
    solver = np.ctypeslib.load_library("libs/euler_steady_solver.so", ".")

    class MeshParam(Structure):
        _fields_ = [("nelem", c_int),
                    ("nnode", c_int),
                    ("niedge", c_int),
                    ("nbedge", c_int),
                    ("nbgroup", c_int)]

    mesh_param = MeshParam(np.shape(mesh.E)[0], np.shape(mesh.V)[0], np.shape(mesh.I2E)[0],
                            np.shape(mesh.B2E)[0], np.shape(mesh.Bname)[0])
    

    # Passing arrays from numpy array to C
    c_E = ((c_int * 3) * np.shape(mesh.E)[0])()
    c_V = ((c_double * 2) * np.shape(mesh.V)[0])()
    c_I2E = ((c_int * 4) * np.shape(mesh.I2E)[0])()
    c_B2E = ((c_int * 3) * np.shape(mesh.B2E)[0])()
    c_In = ((c_double * 2) * np.shape(mesh.In)[0])()
    c_Bn = ((c_double * 2) * np.shape(mesh.Bn)[0])()
    c_Area = (c_double * np.shape(mesh.Area)[0])()
    c_state_vectors = ((c_double * 4) * np.shape(mesh.E)[0])()
    c_gradu = (((c_double * 2) * 4) * np.shape(mesh.E)[0])()

    # Set up the free strem initial states
    for i in range(np.shape(mesh.E)[0]):
        c_state_vectors[i][0] = state_vectors[i][0]
        c_state_vectors[i][1] = state_vectors[i][1]
        c_state_vectors[i][2] = state_vectors[i][2]
        c_state_vectors[i][3] = state_vectors[i][3]

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

    solver.calculate_gradient(byref(mesh_param),
                            c_E, c_V, c_I2E, c_B2E, c_In, c_Bn, c_Area,
                            c_state_vectors, c_gradu)
    
    gradu = np.zeros([np.shape(mesh.E)[0], 4, 2]) 

    for i in range(np.shape(mesh.E)[0]):
        for j in range(4):
            for k in range(2):
                gradu[i, j, k] = c_gradu[i][j][k]

    return gradu


    

    

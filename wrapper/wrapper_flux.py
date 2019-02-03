from ctypes import *
import numpy as np

def c_roe_2d(py_uL, py_uR, py_n, py_gamma):
    flux = np.ctypeslib.load_library("libs/flux.so", ".")
    c_uL = (c_double * 4)(py_uL[0], py_uL[1], py_uL[2], py_uL[3])
    c_uR = (c_double * 4)(py_uR[0], py_uR[1], py_uR[2], py_uR[3])
    c_n = (c_double * 2)(py_n[0], py_n[1])
    c_gamma = c_double(py_gamma)
    c_F_hat = (c_double * 4)()
    c_mws = c_double()
    flux.roe_2d(c_uL, c_uR, c_n, c_gamma, c_F_hat, byref(c_mws))
    arr_F_hat = np.zeros([4, 1])
    for i in range(4):
        arr_F_hat[i] = c_F_hat[i]
    return arr_F_hat, c_mws.value

def c_projection_2d(py_F, py_n):
    flux = np.ctypeslib.load_library("libs/flux.so", ".")
    c_F = ((c_double * 2) * 4)()
    c_n = (c_double * 2)(py_n[0], py_n[1], py_n[2], py_n[3])
    c_F_hat = (c_double * 4)()
    flux.projection_2d(c_F, c_n, c_F_hat)
    F_hat_arr = np.zeros([4, 1])
    for i in range(4):
        F_hat_arr[i] = c_F_hat[i]
    return F_hat_arr

def c_flux_function_2d(py_u, py_gamma):
    flux = np.ctypeslib.load_library("libs/flux.so", ".")
    c_F = ((c_double * 2) * 4)()
    c_u = (c_double * 4)(py_u[0], py_u[1], py_u[2], py_u[3])
    c_gamma = c_double(py_gamma)
    F_arr = np.zeros([4, 2])
    flux.flux_function_2d(c_u, c_gamma, c_F)
    for i in range(4):
        for j in range(2):
            F_arr[i][j] = c_F[i][j]
    return F_arr


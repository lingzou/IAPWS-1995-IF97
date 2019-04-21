import ctypes
from ctypes import cdll
from ctypes import *
import numpy as np
import matplotlib.pyplot as plt
import os


#lib_IF97 = ctypes.CDLL('lib_IF97.so')
lib_IF97 = cdll.LoadLibrary('../libIF97.so')
lib_IF97.B23_p_from_T.argtypes = [c_double]
lib_IF97.B23_p_from_T.restype = c_double

lib_IF97.p_sat_from_T.argtypes = [c_double]
lib_IF97.p_sat_from_T.restype = c_double
lib_IF97.T_sat_from_p.argtypes = [c_double]
lib_IF97.T_sat_from_p.restype = c_double

'''
rho_sat_data = np.loadtxt(os.path.join(os.pardir, 'tables', 'rho_sat_from_T.dat'))
T_array = rho_sat_data[:, 0]
rho_l_sat_array = rho_sat_data[:, 1]
rho_g_sat_array = rho_sat_data[:, 2]

fig, ax = plt.subplots()
ax.plot(T_array, rho_l_sat_array, color = "b", marker = ".", ls="-")
ax.plot(T_array, rho_g_sat_array, color = "k", marker = ".", ls="-")
plt.show()'''

print lib_IF97.B23_p_from_T(0.62315e3)
print lib_IF97.p_sat_from_T(300.0)
print lib_IF97.T_sat_from_p(10.0e6)

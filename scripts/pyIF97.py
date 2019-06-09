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
lib_IF97.R3_rho_l_sat_from_T_ITER.argtypes = [c_double]
lib_IF97.R3_rho_l_sat_from_T_ITER.restype = c_double

sat_data = np.loadtxt(os.path.join(os.pardir, 'tables', 'sat_line.dat'), skiprows = 1)
T_array = sat_data[:, 0]
p_sat_array = sat_data[:, 1]
#rho_g_sat_array = sat_data[:, 2]
T23_array = np.arange(623.15, 863.15, 1.0)
p23_array = np.arange(623.15, 863.15, 1.0)
for i in xrange(0, T23_array.shape[0]) :
  p23_array[i] = lib_IF97.B23_p_from_T(T23_array[i])

fig, ax = plt.subplots()
ax.plot([273.15, 273.15], [100.0, 100.0e6], color = "r", marker = "None", ls = "-")
ax.plot([273.15, 1073.15], [100.0e6, 100.0e6], color = "r", marker = "None", ls = "-")
ax.plot([1073.15, 1073.15], [100.0e6, 100.0], color = "r", marker = "None", ls = "-")
ax.plot([623.15, 623.15], [lib_IF97.p_sat_from_T(623.15), 100.0e6], color = "k", marker = "None", ls = "-")
ax.plot(T23_array, p23_array, color = "k", marker = "None", ls="-")
ax.plot(T_array, p_sat_array, color = "b", marker = "None", ls="-")
#ax.set_yscale('log')
plt.show()

print lib_IF97.B23_p_from_T(0.62315e3)
print lib_IF97.p_sat_from_T(300.0)
print lib_IF97.T_sat_from_p(10.0e6)

print lib_IF97.R3_rho_l_sat_from_T_ITER(623.149)
print lib_IF97.R3_rho_l_sat_from_T_ITER(623.15)
print lib_IF97.R3_rho_l_sat_from_T_ITER(623.151)

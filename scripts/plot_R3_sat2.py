import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.tri as mtri
from scipy.interpolate import CubicSpline
import math


data = np.loadtxt('R3_rho_sat.dat')
fig, ax = plt.subplots()

T_array = data[:, 0]
rho_l_sat_array = data[:, 1]
rho_g_sat_array = data[:, 2]

cs_l = CubicSpline(T_array, rho_l_sat_array)
cs_g = CubicSpline(T_array, rho_g_sat_array)

TT = np.arange(T_array[40], T_array[54], 0.0001)
TT_last = np.arange(51, dtype=np.float64)
rho_last = np.arange(51, dtype=np.float64)
Tc = T_array[54]
T2 = T_array[53]
dT = (Tc-T2)/50
rho_c = rho_l_sat_array[54]
rho_2 = rho_l_sat_array[53]
print rho_2
#a = pow(Tc - T2, 0.5) / (rho_2 - rho_c)
a = (T2 - Tc) / (rho_2 - rho_c) / (rho_2 - rho_c)

print TT_last
print rho_last
print 'a=', a

for i in xrange(0, 50) :
  TT_last[i] = 647.0955 + dT * i
  print TT_last[i]
  rho_last[i] = rho_c + pow((TT_last[i] - Tc) / a, 0.5)
  print rho_last[i]
TT_last[50] = 647.096
rho_last[50] = 322.0



ax.plot(T_array, rho_l_sat_array, color = "#DC143C", marker = ".", ls="")
ax.plot(T_array, rho_g_sat_array, color = "k", marker = ".", ls="")

ax.plot(TT, cs_l(TT), color = "#DC143C", marker = "None", ls="-")
ax.plot(TT, cs_g(TT), color = "k", marker = "None", ls="-")

ax.plot(TT_last, rho_last, color = "k", marker = "+", ls="-")

plt.show()

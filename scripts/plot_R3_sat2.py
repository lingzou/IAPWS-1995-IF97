import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
import math
import os

#data = np.loadtxt(os.path.join(os.pardir, 'UnitTest', 'R3_rho_sat.dat'))
data = np.loadtxt(os.path.join(os.pardir, 'R3_rho_sat.dat'))
fig, ax = plt.subplots()

T_array = data[:, 0]
rho_l_sat_array = data[:, 1]
rho_g_sat_array = data[:, 2]

T_array_m1 = T_array[:-1]
rho_l_sat_array_m1 = rho_l_sat_array[:-1]
rho_g_sat_array_m1 = rho_g_sat_array[:-1]

#cs_l = CubicSpline(T_array, rho_l_sat_array)
#cs_g = CubicSpline(T_array, rho_g_sat_array)
cs_l = CubicSpline(T_array_m1, rho_l_sat_array_m1)
cs_g = CubicSpline(T_array_m1, rho_g_sat_array_m1)

TT = np.arange(T_array[40], T_array[54], 0.0001)
TT_last = np.arange(51, dtype=np.float64)
rho_l_last = np.arange(51, dtype=np.float64)
rho_g_last = np.arange(51, dtype=np.float64)
Tc = T_array[54]
T2 = T_array[53]
dT = (Tc-T2)/50
rho_c = rho_l_sat_array[54]
rho_2 = rho_l_sat_array[53]
rho_1 = rho_g_sat_array[53]
print rho_2, rho_1
#a = pow(Tc - T2, 0.5) / (rho_2 - rho_c)
a2 = (T2 - Tc) / (rho_2 - rho_c) / (rho_2 - rho_c)
a1 = (T2 - Tc) / (rho_1 - rho_c) / (rho_1 - rho_c)

print TT_last
print rho_l_last
print rho_g_last
print 'a2=', a2
print 'a1=', a1

for i in xrange(0, 50) :
  TT_last[i] = 647.0955 + dT * i
  print TT_last[i]
  rho_l_last[i] = rho_c + pow((TT_last[i] - Tc) / a2, 0.5)
  print rho_l_last[i]
  rho_g_last[i] = rho_c - pow((TT_last[i] - Tc) / a1, 0.5)
  print rho_g_last[i]
TT_last[50] = 647.096
rho_l_last[50] = 322.0
rho_g_last[50] = 322.0



ax.plot(T_array, rho_l_sat_array, color = "#DC143C", marker = ".", ls="-")
ax.plot(T_array, rho_g_sat_array, color = "k", marker = ".", ls="-")

ax.plot(TT, cs_l(TT), color = "#DC143C", marker = "None", ls="-")
ax.plot(TT, cs_g(TT), color = "k", marker = "None", ls="-")

#ax.plot(TT_last, cs_l(TT_last), color = "#DC143C", marker = "None", ls="-")
#ax.plot(TT_last, cs_g(TT_last), color = "k", marker = "None", ls="-")

ax.plot(TT_last, rho_l_last, color = "k", marker = "+", ls="-")
ax.plot(TT_last, rho_g_last, color = "b", marker = "+", ls="-")

plt.show()

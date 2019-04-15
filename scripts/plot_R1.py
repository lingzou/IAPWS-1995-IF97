import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.tri as mtri

data = np.loadtxt('R1_rho_p_T.dat')

T_array = data[:,0]
p_array = data[:,1]
rho_array = data[:,2]

fig = plt.figure()
ax = Axes3D(fig)
#ax = plt.axes(projection='3d')
ax.scatter3D(T_array, p_array, rho_array, c=rho_array, cmap='Greens');
plt.show()

#fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
#
#ax.scatter(T_array, p_array, marker=".", c="#DC143C", edgecolors="black", s=100)
#ax.set_xlabel('X')
#ax.set_ylabel('Y')
#
#plt.show()
#
#triang = mtri.Triangulation(T_array, p_array)
#
#fig2 = plt.figure()
#ax = fig2.add_subplot(1,1,1)
#
#ax.triplot(triang, c="#D3D3D3", marker='.', markerfacecolor="#DC143C", markeredgecolor="black", markersize=1)
#
#ax.set_xlabel('X')
#ax.set_ylabel('Y')
#plt.show()

#print data

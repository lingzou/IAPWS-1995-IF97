import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.tri as mtri

for i in xrange(0, 51) :
  data = np.loadtxt('../output/' + str(i) + '.dat')

  rho1_array = data[0:51,0]
  val1_array = data[0:51,1]

  rho2_array = data[51:102,0]
  val2_array = data[51:102,1]

  rho3_array = data[102:153,0]
  val3_array = data[102:153,1]

  rho_sat_array = data[153:155,0]
  val_sat_array = data[153:155,1]

  fig, ax = plt.subplots()
  #ax = fig.add_subplot()

  #plt.plot(rho_array, val_array, marker=".", c="#DC143C", edgecolors="black")
  ax.plot(rho1_array, val1_array, color = "#DC143C", marker = None, ls="-")
  #ax.plot(rho2_array, val2_array, color = 'k', marker = ".", ls="")
  ax.plot(rho2_array, val2_array, color = "#DC143C", marker = None, ls="-")
  ax.plot(rho3_array, val3_array, color = "#DC143C", marker = None, ls="-")
  ax.plot(rho_sat_array, val_sat_array, color = "b", marker = "x", ls="")
  #ax.plot([data[0,0], data[152,0]], [0.0, 0.0], color = "k", marker = None, ls="-", lw=0.5)
  ax.plot([0, 1000], [0, 0], color = "k", marker = None, ls="-", lw=0.5)
  ax.set_xlabel(r'$\rho$ [kg/m$^3$]')
  ax.set_ylabel(r'$\frac{p_s}{RT\rho}-\delta \phi_\delta(\delta, \tau)$')
  ax.set_xlim([data[0,0], data[152,0]])
  plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

  #ax.autoscale(enable = True, axis = "y", tight = False)

  plt.show()
  #plt.draw()
  #plt.savefig('three-roots.pdf', dpi=600)

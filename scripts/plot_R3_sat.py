import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.tri as mtri

for i in xrange(0, 54) :
  data = np.loadtxt('output/' + str(i) + '.dat')

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
  ax.plot(rho1_array, val1_array, color = "#DC143C", marker = ".", ls="")
  ax.plot(rho2_array, val2_array, color = 'k', marker = ".", ls="")
  ax.plot(rho3_array, val3_array, color = "#DC143C", marker = ".", ls="")
  ax.plot(rho_sat_array, val_sat_array, color = "y", marker = "+", ls="-")
  ax.plot([data[0,0], data[152,0]], [0.0, 0.0], color = "k", marker = ".", ls="-")
  #plt.set_xlabel('Rho')
  #plt.set_ylabel('Val')

  ax.autoscale(enable = True, axis = "y", tight = False)

  plt.show()

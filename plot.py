import os, sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import scipy.interpolate as interp
from netCDF4 import Dataset
import operator

in_file =  str(sys.argv[1])
var1 =  str(sys.argv[2])
var2 =  str(sys.argv[3])
idx = int(sys.argv[4])
ncfile = Dataset(in_file, 'r', format='NETCDF3')

#Read relevant data from TRANSP file
x = np.array(ncfile.variables[var1][:])
y = np.array(ncfile.variables[var2][:])

if var1 == 'TIME' or var1 == 'TIME3':
  plt.plot(x[:,idx],y[:,idx])
  plt.show()
else:
  plt.plot(x[idx,:],y[idx,:])
  plt.show()
  
